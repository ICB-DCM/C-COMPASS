"""Create synthetic datasets for C-COMPASS testing."""

import re
from collections import Counter

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from ccompass.core import (
    IDENTIFIER,
    KEEP,
    NA,
)


class SyntheticDataConfig(BaseModel):
    """Configuration for synthetic data generation."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: number of simulated replicates per condition
    replicates: int = 3
    #: number of simulated conditions
    conditions: int = 2
    #: number of simulated compartments
    num_compartments: int = 5
    #: number of fractions per gradient
    fractions: int = 14

    #: peak width in fractions (range)
    spread: list[int] = [3, 4]
    #: peak height as intensities (range)
    intensities: list[int] = [5000, 8000]
    #: increase for more intensity-variance across a single peaks
    intensity_variance: int = 50
    #: increase this value for more accurate peaks
    peak_accuracy: int = 200
    #: number of marker proteins per compartment (range)
    markers: list[int] = [50, 100]
    #: number of proteins per compartment with one unknown localization (range)
    unknown_single: list[int] = [50, 100]
    #: number of proteins per compartment with a random second localization (range)
    unknown_double: list[int] = [10, 20]
    #: possible mixing ratios for double localizations
    ratios_double: list[tuple[int, int]] = [(75, 25), (50, 50)]
    #: number of proteins per compartment with two random additional localizations (range)
    unknown_triple: list[int] = [3, 8]
    #: possible mixing ratios for triple localizations
    ratios_triple: list[tuple[int, int, int]] = [(50, 25, 25)]
    #: chance for a protein to be missing in single replicates
    missing_rep: float = 0.04
    #: chance for a protein to be missing in single conditions
    missing_cond: float = 0.02
    #: probability for relocalization
    reloc_rate: float = 0.1

    #: protein ID column name
    protein_id_col: str = "ProteinName"
    #: gene ID column name
    gene_id_col: str = "GeneName"
    #: class ID column name
    class_id_col: str = "Marker"

    #: random number generator seed
    rng_seed: int = 1
    #: random number generator
    rng: np.random.Generator = Field(
        default_factory=lambda data: np.random.default_rng(data["rng_seed"])
    )


def create_compspecs(c: SyntheticDataConfig) -> dict:
    """Create compartment specifications."""
    comp_specs = {}
    rng = c.rng

    for i in range(c.num_compartments):
        name = f"Compartment{i + 1}"
        comp_specs[name] = {}

        bins = np.linspace(1, c.fractions, c.num_compartments + 1)
        middle = rng.integers(int(bins[i]), int(bins[i + 1]))
        comp_specs[name]["middle"] = middle

        sigma = rng.uniform(c.spread[0] / 2, c.spread[1] / 2)
        comp_specs[name]["sigma"] = sigma

        height = rng.uniform(c.intensities[0], c.intensities[1])
        comp_specs[name]["height"] = height

        comp_specs[name]["number_marker"] = int(
            rng.uniform(c.markers[0], c.markers[1])
        )
        comp_specs[name]["number_single"] = int(
            rng.uniform(c.unknown_single[0], c.unknown_single[1])
        )
        comp_specs[name]["number_double"] = int(
            rng.uniform(c.unknown_double[0], c.unknown_double[1])
        )
        comp_specs[name]["number_triple"] = int(
            rng.uniform(c.unknown_triple[0], c.unknown_triple[1])
        )
    return comp_specs


def reflect_distribution(distribution, lower, upper):
    reflected = np.copy(distribution)
    reflected[distribution < lower] = lower
    reflected[distribution > upper] = upper
    return reflected


def create_profile(middle, sigma, height, c: SyntheticDataConfig):
    """Create a single profile."""
    rng = c.rng
    num_samples = int(rng.normal(c.peak_accuracy, c.peak_accuracy / 10))
    random_values = rng.normal(middle, sigma, num_samples)

    reflected_values = reflect_distribution(random_values, 1, c.fractions)
    discrete_values = np.clip(
        np.round(reflected_values), 1, c.fractions
    ).astype(int)
    value_counts = Counter(discrete_values)
    value_counts_x = np.array(
        [
            value_counts[i] if i in value_counts else 0
            for i in range(1, c.fractions + 1)
        ]
    )

    factor = height / max(value_counts_x)
    value_counts_scaled = np.copy(value_counts_x)
    for j in range(c.fractions):
        factor_rand = float(rng.normal(factor, height / c.intensity_variance))
        new_value = value_counts_x[j] * factor_rand
        if new_value >= 0:
            value_counts_scaled[j] = new_value
        else:
            while new_value < 0:
                factor_rand = float(rng.normal(factor, height / 45))
                new_value = value_counts_x[j] * factor_rand
            value_counts_scaled[j] = new_value
    return value_counts_scaled


def create_profiles(c: SyntheticDataConfig):
    rng = c.rng
    comp_specs = {}
    for cond in range(c.conditions):
        comp_specs[cond] = create_compspecs(c)
    comp_list = [f"Compartment{i}" for i in range(1, c.num_compartments + 1)]
    data_columns = [c.protein_id_col, c.gene_id_col]
    for cond in range(c.conditions):
        for rep in range(c.replicates):
            for fract in range(c.fractions):
                data_columns.append(
                    f"Con{cond + 1}_Rep{rep + 1}_Fr{str(fract + 1).zfill(2)}"
                )
        for comp in comp_specs[cond]:
            data_columns.append(f"Amount_{cond}_{comp}")
    data_columns.append(c.class_id_col)

    # CREATE MARKER PROFILES:
    count = 0
    all_profiles = []
    for comp in comp_list:
        for m in range(comp_specs[cond][comp]["number_marker"]):
            count = count + 1
            prot_name = f"Prot{count}"
            gene_name = f"Gene{count}"
            profile_conc = [prot_name, gene_name]
            for cond in range(c.conditions):
                empty_condition = rng.random() < c.missing_cond
                if empty_condition:
                    profile_conc.extend(
                        (c.fractions * c.replicates + len(comp_list))
                        * [np.nan]
                    )
                else:
                    specs = comp_specs[cond]
                    location = [0 if key != comp else 1 for key in specs]
                    middle_1 = specs[comp]["middle"]
                    sigma_1 = specs[comp]["sigma"]
                    height_1 = specs[comp]["height"]
                    for rep in range(c.replicates):
                        empty_replicate = rng.random() < c.missing_rep
                        if empty_replicate:
                            profile_conc.extend(c.fractions * [np.nan])
                        else:
                            profile_1 = create_profile(
                                middle_1, sigma_1, height_1, c
                            )
                            profile_conc.extend(profile_1.astype(float))
                    profile_conc.extend(location)
            profile_conc.append(comp)
            all_profiles.append(profile_conc)

    # CREATE UNKNOWN SINGLE LOCALIZATIONS:
    for comp in comp_list:
        comp_others = [c for c in comp_list if c != comp]
        for s in range(comp_specs[cond][comp]["number_single"]):
            count = count + 1
            prot_name = f"Prot{count}"
            gene_name = f"Gene{count}"
            profile_conc = [prot_name, gene_name]
            for cond in range(c.conditions):
                empty_condition = rng.random() < c.missing_cond
                if empty_condition:
                    profile_conc.extend(
                        (c.fractions * c.replicates + len(comp_list))
                        * [np.nan]
                    )
                else:
                    specs = comp_specs[cond]
                    reloc = rng.random() < c.reloc_rate
                    if reloc:
                        cc = rng.choice(comp_others)
                    else:
                        cc = comp
                    location = [0 if key != cc else 1 for key in specs]
                    middle_1 = specs[cc]["middle"]
                    sigma_1 = specs[cc]["sigma"]
                    height_1 = specs[cc]["height"]
                    for rep in range(c.replicates):
                        empty_replicate = rng.random() < c.missing_rep
                        if empty_replicate:
                            profile_conc.extend(c.fractions * [np.nan])
                        else:
                            profile_1 = create_profile(
                                middle_1, sigma_1, height_1, c
                            )
                            profile_conc.extend(profile_1.astype(float))
                    profile_conc.extend(location)
            profile_conc.append(np.nan)
            all_profiles.append(profile_conc)

    # CREATE UNKNOWN DOUBLE LOCALIZATIONS:
    for comp in comp_list:
        comp_others = [c for c in comp_list if c != comp]
        assert (
            c.num_compartments >= 2
            or comp_specs[cond][comp]["number_double"] == 0
        )
        for d in range(comp_specs[cond][comp]["number_double"]):
            count = count + 1
            prot_name = f"Prot{count}"
            gene_name = f"Gene{count}"
            profile_conc = [prot_name, gene_name]
            for cond in range(c.conditions):
                empty_condition = rng.random() < c.missing_cond
                if empty_condition:
                    profile_conc.extend(
                        (c.fractions * c.replicates + len(comp_list))
                        * [np.nan]
                    )
                else:
                    specs = comp_specs[cond]
                    reloc = rng.random() < c.reloc_rate
                    if reloc:
                        c_1 = rng.choice(comp_others)
                        c_others = [co for co in comp_list if co != c_1]
                        c_2 = rng.choice(c_others)
                    else:
                        c_1 = comp
                        c_2 = rng.choice(comp_others)

                    middle_1 = specs[c_1]["middle"]
                    sigma_1 = specs[c_1]["sigma"]
                    height_1 = specs[c_1]["height"]
                    middle_2 = specs[c_2]["middle"]
                    sigma_2 = specs[c_2]["sigma"]
                    height_2 = specs[c_2]["height"]

                    ratio = rng.choice(c.ratios_double)
                    location = [
                        0
                        if key != c_1 and key != c_2
                        else ratio[0] / 100
                        if key == c_1
                        else ratio[1] / 100
                        for key in specs
                    ]

                    for rep in range(c.replicates):
                        empty_replicate = rng.random() < c.missing_rep
                        if empty_replicate:
                            profile_conc.extend(c.fractions * [np.nan])
                        else:
                            profile_1 = create_profile(
                                middle_1, sigma_1, height_1, c
                            ).astype(float)
                            profile_2 = create_profile(
                                middle_2, sigma_2, height_2, c
                            ).astype(float)
                            profile_combined = (ratio[0] / 100 * profile_1) + (
                                ratio[1] / 100 * profile_2
                            )
                            profile_conc.extend(profile_combined.astype(float))
                    profile_conc.extend(location)
            profile_conc.append(np.nan)
            all_profiles.append(profile_conc)

    # CREATE UNKNOWN TRIPLE LOCALIZATIONS:
    assert (
        c.num_compartments >= 3 or comp_specs[cond][comp]["number_triple"] == 0
    )
    for comp in comp_list:
        comp_others = [c for c in comp_list if c != comp]
        for d in range(comp_specs[cond][comp]["number_triple"]):
            count = count + 1
            prot_name = f"Prot{count}"
            gene_name = f"Gene{count}"
            profile_conc = [prot_name, gene_name]
            for cond in range(c.conditions):
                empty_condition = rng.random() < c.missing_cond
                if empty_condition:
                    profile_conc.extend(
                        (c.fractions * c.replicates + len(comp_list))
                        * [np.nan]
                    )
                else:
                    specs = comp_specs[cond]
                    reloc = rng.random() < c.reloc_rate
                    if reloc:
                        c_1 = rng.choice(comp_others)
                        c_others = [co for co in comp_list if co != c_1]
                        c_2 = rng.choice(c_others)
                        c_others = [
                            co for co in comp_list if co != c_1 and co != c_2
                        ]
                        c_3 = rng.choice(c_others)
                    else:
                        c_1 = comp
                        c_2 = rng.choice(comp_others)
                        c_others = [co for co in comp_others if co != c_2]
                        c_3 = rng.choice(c_others)

                    middle_1 = specs[c_1]["middle"]
                    sigma_1 = specs[c_1]["sigma"]
                    height_1 = specs[c_1]["height"]
                    middle_2 = specs[c_2]["middle"]
                    sigma_2 = specs[c_2]["sigma"]
                    height_2 = specs[c_2]["height"]
                    middle_3 = specs[c_3]["middle"]
                    sigma_3 = specs[c_3]["sigma"]
                    height_3 = specs[c_3]["height"]

                    ratio = rng.choice(c.ratios_triple)
                    location = [
                        0
                        if key != c_1 and key != c_2 and key != c_3
                        else ratio[0] / 100
                        if key == c_1
                        else ratio[1] / 100
                        if key == c_2
                        else ratio[2] / 100
                        for key in specs
                    ]

                    for rep in range(c.replicates):
                        empty_replicate = rng.random() < c.missing_rep
                        if empty_replicate:
                            profile_conc.extend(c.fractions * [np.nan])
                        else:
                            profile_1 = create_profile(
                                middle_1, sigma_1, height_1, c
                            ).astype(float)
                            profile_2 = create_profile(
                                middle_2, sigma_2, height_2, c
                            ).astype(float)
                            profile_3 = create_profile(
                                middle_3, sigma_3, height_3, c
                            ).astype(float)
                            profile_combined = (
                                (ratio[0] / 100 * profile_1)
                                + (ratio[1] / 100 * profile_2)
                                + (ratio[2] / 100 * profile_3)
                            )
                            profile_conc.extend(profile_combined.astype(float))
                    profile_conc.extend(location)
            profile_conc.append(np.nan)
            all_profiles.append(profile_conc)

    dataset = pd.DataFrame(all_profiles, columns=data_columns)
    dataset_shuffled = dataset.sample(frac=1, random_state=rng).reset_index(
        drop=True
    )

    markerset = dataset[[c.gene_id_col, c.class_id_col]].dropna(
        subset=[c.class_id_col]
    )

    return dataset_shuffled, markerset


def total_proteome(
    proteins: list[str], c: SyntheticDataConfig
) -> pd.DataFrame:
    """Generate total proteome data."""
    tp_intensities = [7000, 10000]
    variance = 500
    regulated = 20  # precentage of proteins with changing expression level
    changes = [0.1, 10]  # possible fold changes for regulated proteins (range)
    rng = c.rng

    tp_columns = [c.protein_id_col]
    for cond in range(c.conditions):
        for rep in range(c.replicates):
            colname = f"Con{cond + 1}_Rep{rep + 1}"
            tp_columns.append(colname)
        tp_columns.append(f"RelativeRegulation_Con{cond + 1}")

    all_values = []
    for prot in proteins:
        height = rng.uniform(tp_intensities[0], tp_intensities[1])
        values = [prot]
        changing = rng.random() < regulated / 100

        if changing:
            fc = rng.uniform(changes[0], changes[1])
        else:
            fc = rng.normal(1, 0.2)

        for cond in range(c.conditions):
            height_cond = height * fc
            variance_cond = variance * fc
            for rep in range(c.replicates):
                values.append(rng.normal(height_cond, variance_cond))
            values.append(fc)
        all_values.append(values)

    return pd.DataFrame(all_values, columns=tp_columns)


# regexes to parse column IDs
fract_id_rx = re.compile(
    r"(?P<condition>Con\d+)_Rep(?P<replicate>\d+)_Fr(?P<fraction>\d+)"
)
tp_id_rx = re.compile(r"(?P<condition>Con\d+)_Rep(?P<replicate>\d+)")


def fract_col_id_to_row(col_id: str, c: SyntheticDataConfig) -> list:
    """Convert fractionation data column id to fractionation table rows."""
    if col_id == c.protein_id_col:
        return [col_id, IDENTIFIER, NA, NA]
    if col_id == c.gene_id_col:
        return [col_id, KEEP, NA, NA]

    if not (match := fract_id_rx.match(col_id)):
        raise ValueError(f"Invalid fractionation ID: {col_id}")

    condition = match["condition"]
    replicate = int(match["replicate"])
    fraction = int(match["fraction"])
    return [col_id, condition, replicate, fraction]


def tp_col_id_to_row(col_id: str, c: SyntheticDataConfig) -> list:
    """Convert total proteome data column id to total proteome table rows."""
    if col_id == c.protein_id_col:
        return [col_id, IDENTIFIER]

    if not (match := tp_id_rx.match(col_id)):
        raise ValueError(f"Invalid total proteome ID: {col_id}")

    condition = match["condition"]
    return [col_id, condition]


def main():
    # filenames for the synthetic data
    filename_fract = "sim_Fractionation.txt"
    filename_marker = "sim_Markerlist.txt"
    filename_tp = "sim_TotalProteome.txt"

    c = SyntheticDataConfig()
    dataset, markerset = create_profiles(c)
    dataset_tp = total_proteome(proteins=list(dataset[c.protein_id_col]), c=c)

    dataset.to_csv(filename_fract, sep="\t", index=False)
    markerset.to_csv(filename_marker, sep="\t", index=False)
    dataset_tp.to_csv(filename_tp, sep="\t", index=False)


if __name__ == "__main__":
    main()
