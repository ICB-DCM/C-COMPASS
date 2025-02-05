"""Core classes and functions for the ccompass package."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, field_serializer

from . import config_filepath


class AppSettings(BaseModel):
    """Settings for the C-COMPASS application"""

    #: The directory that was last used to load/save a session
    last_session_dir: Path = Path.home()

    @field_serializer("last_session_dir")
    def serialize_last_session_dir(self, value: Path) -> str:
        return str(value)

    @classmethod
    def load(cls, filepath: Path = None):
        """Load the settings from a file."""
        import yaml

        if filepath is None:
            filepath = config_filepath

        if not filepath.exists():
            return cls()

        with open(filepath) as f:
            data = yaml.safe_load(f) or {}
            return cls(**data)

    def save(self, filepath: Path = None):
        """Save the settings to a file."""
        import yaml

        if filepath is None:
            filepath = config_filepath

        filepath.parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, "w") as f:
            yaml.safe_dump(self.model_dump(), f)


class NeuralNetworkParametersModel(BaseModel):
    """Hyperparameters for the neural network."""

    #: Perform upsampling?
    upsampling: bool = True
    #: Method for upsampling
    upsampling_method: Literal[
        "none", "noised", "average", "noisedaverage"
    ] = "noisedaverage"
    #: Noise level for upsampling (standard deviations)
    upsampling_noise: float = 2
    #: Auto-encoder type
    AE: Literal["none", "lite", "full", "full_lite"] = "none"
    # FIXME: unused
    AE_activation: Literal["relu", "leakyrelu"] = "leakyrelu"
    # FIXME: unused
    AE_out: Literal["sigmoid", "relu", "softmax", "leakyrelu"] = "sigmoid"
    # FIXME: unused
    AE_epochs: int = 20
    #: Perform SVM filtering?
    svm_filter: bool = False
    #: ...
    # FIXME: can be "none"; == 0?!
    mixed_part: int | str = 4
    #: The fraction of the mixed batch to use (0-1)
    mixed_batch: float = 0.05
    #: Long or short optimization?
    NN_optimization: Literal["short", "long"] = "long"
    #: Neural network activation function
    NN_activation: Literal["relu", "leakyrelu"] = "relu"
    #: Neural network class layer activation function
    class_activation: Literal["sigmoid", "softmax", "linear"] = "linear"
    #: Neural network training loss function
    class_loss: Literal["binary_crossentropy", "mean_squared_error"] = (
        "mean_squared_error"
    )
    #: FIXME: unused
    regularization: Literal["none", "l1", "l2", "elastic"] = "none"
    #: Optimizers to include in the hyperparameter search
    optimizers: list[Literal["adam", "rmsprop", "sgd"]] = [
        "adam",
        "rmsprop",
        "sgd",
    ]
    #: Number of epochs for the neural network training
    NN_epochs: int = 20
    #: ...
    rounds: int = 3
    #: ...
    subrounds: int = 10
    #: Percentile threshold for ... ?
    reliability: int = 95


class SessionStatusModel(BaseModel):
    """Keeps track of the different analysis steps that have been completed."""

    fractionation_data: bool = False
    tp_data: bool = False
    lipidome_data: bool = False
    lipidome_total: bool = False
    marker_file: bool = False
    marker_matched: bool = False
    training: bool = False
    proteome_prediction: bool = False
    lipidome_prediction: bool = False
    comparison_global: bool = False
    comparison_class: bool = False


def fract_default():
    """Default settings for fractionation data processing."""
    params_default = {
        "class": {
            "scale1": [
                True,
                "area",
            ],
            "corrfilter": False,
            "scale2": [False, "area"],
            "zeros": True,
            "combination": "separate",
        },
        "vis": {
            "scale1": [
                True,
                "minmax",
            ],
            "corrfilter": False,
            "scale2": [True, "minmax"],
            "zeros": True,
            "combination": "median",
        },
        "global": {
            "missing": [True, "1"],
            "minrep": [True, "2"],
            "outcorr": False,
        },
    }
    return params_default


# type annotations

# A condition ID
ConditionId = str
# Path to a file
Filepath = str
# Condition + replicate ID: "{condition}_Rep.{replicate}"
ConditionReplicate = str


class SessionModel(BaseModel):
    """Data for a C-COMPASS session."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: Filepaths for fractionation data
    fract_paths: list[Filepath] = []
    #: Fractionation column assignments
    #  filepath => [column ID, condition, replicate, fraction]
    fract_tables: dict[Filepath, list[list[int | str]]] = {}
    #: ??
    fract_pos: dict[Filepath, list[int]] = {}
    #: Fractionation input files: filepath => DataFrame
    fract_indata: dict[Filepath, pd.DataFrame] = {}
    #: Fractionation data for classification and visualization
    #  One DataFrame for each condition x replicate
    #  ("{condition}_Rep.{replicate}")
    fract_data: dict[ConditionReplicate, dict[str, pd.DataFrame]] = {
        "class": {},
        "vis": {},
    }
    #: ??
    #  for visualization and classification, each containing one DataFrame
    #  per condition with columns "{condition}_std_Fr.{fraction}"
    fract_std: dict[
        Literal["class", "vis"], dict[ConditionId, pd.DataFrame]
    ] = {"class": {}, "vis": {}}
    #: ??
    #  *something* => "{condition}" => "Rep.{replicate}" => DataFrame
    fract_intermediate: dict[str, dict[str, dict[str, pd.DataFrame]]] = {}
    #: Identifier column of each fractionation data table:
    #   filepath => column id
    fract_identifiers: dict[Filepath, str] = {}
    #: Addition ("keep") columns from the fractionation data
    #  column ID => DataFrame
    fract_info: dict[str, pd.DataFrame] = {}
    #: Fractionation preprocessing parameters.
    #  global/classification/visualization
    #  "global"|"class"|"vis" => option => value
    fract_preparams: dict[str, dict[str, Any]] = fract_default()
    #: Conditions in the fractionation data, including "[KEEP]"
    fract_conditions: list[str] = []
    #: Fractionation data for the different conditions x replicates
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_full: dict[ConditionReplicate, pd.DataFrame] = {}
    #: Fractionation data after upsampling
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_full_up: dict[ConditionReplicate, pd.DataFrame] = {}
    #: Marker abundance in the different fractions
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_marker: dict[ConditionReplicate, pd.DataFrame] = {}
    #: Marker abundance in the different fractions for visualization
    #  "{condition}_median" => DataFrame
    fract_marker_vis: dict[str, pd.DataFrame] = {}
    #: Marker abundance in the different fractions after upsampling
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_marker_up: dict[ConditionReplicate, pd.DataFrame] = {}
    #: ??
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_mixed_up: dict[ConditionReplicate, pd.DataFrame] = {}
    #: ??
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_test: dict[ConditionReplicate, pd.DataFrame] = {}

    #: Filepaths for total proteome data
    tp_paths: list[Filepath] = []
    #: Total proteome column assignments
    #  filepath => [column ID, condition]
    tp_tables: dict[Filepath, list[list[str]]] = {}
    #: ??
    tp_pos: dict[Filepath, list[int]] = {}
    #: Total proteome input files: filepath => DataFrame
    tp_indata: dict[Filepath, pd.DataFrame] = {}
    #: Total proteome data for the different conditions
    #  One DataFrame for each condition containing all replicates
    #  (column names are "{condition}_Rep.{replicate}")
    tp_data: dict[ConditionReplicate, pd.DataFrame] = {}
    #: ??
    #  *something* => "{condition}" => DataFrame
    tp_intermediate: dict[str, dict[str, pd.DataFrame]] = {}
    #: Identifier column for the total proteome: filepath => column id
    tp_identifiers: dict[Filepath, str] = {}
    #: ??
    tp_icorr: dict = {}
    #: ??
    tp_conditions: list = []
    #: ??
    tp_info: pd.DataFrame = pd.DataFrame()
    #: Total proteome preprocessing parameters
    tp_preparams: dict[str, Any] = {"minrep": 2, "imputation": "normal"}

    #: The user-provided marker files, classes and annotations
    #  filepath => {'table'->pd.DataFrame,
    #  'identifier_col'-> column ID ("key column" in GUI),
    #  'class_col': column ID with class names in the marker file,
    #  'classes': list[str] class names
    #  }
    marker_sets: dict[Filepath, dict[str, Any]] = {}
    #: Options for merging marker sets
    marker_params: dict[str, str] = {"how": "exclude", "what": "unite"}
    #: Mapping of compartment names to class names
    #  nan-values indicate that the compartment is not to be used
    marker_conv: dict[str, str | float] = {}
    #: The consolidated marker list, after merging `marker_sets`
    #  according to `marker_params`, and accounting for renaming
    #  and filtering through `marker_conv`.
    #  "name" (gene name, index) => "class" (class name)
    marker_list: pd.DataFrame = pd.DataFrame()

    #: The column ID of the fractionation DataFrame that is
    #  to be used for matching the markers (`marker_list["name"])
    marker_fractkey: str = "[IDENTIFIER]"

    #: SVM marker prediction
    # "{condition}_Rep.{replicate}" => DataFrame
    #  columns are the fractions + ["class", "svm_prediction", "svm_probability"]
    svm_marker: dict[str, pd.DataFrame] = {}
    #: SVM test data
    # "{condition}_Rep.{replicate}" => DataFrame
    #  columns are the fractions + ["class", "svm_prediction", "svm_probability"]
    svm_test: dict[str, pd.DataFrame] = {}
    #: SVM classification metrics for each condition x replicate
    # "{condition}_Rep.{replicate}" => dict(
    #   "accuracy" -> float,
    #   "precision" -> float,
    #   "recall" -> float,
    #   "f1" -> float,
    #   "confusion" -> pd.DataFrame,
    # )
    svm_metrics: dict[str, dict[str, Any]] = {}

    #: Neural network data
    # "{condition}_Rep.{replicate}" => dict(
    #  {w,W,x,X,y,Y,z,Z}_... => ...
    # )
    learning_xyz: dict[ConditionReplicate, dict[str, Any]] = {}
    #: Nerural network hyperparameters
    NN_params: NeuralNetworkParametersModel = NeuralNetworkParametersModel()

    #: SVM results (?)
    # "{condition}" => dict(
    #  "metrics" -> DataFrame,
    #  "SVM" -> dict("winner_combined" => DataFrame,
    #                "prob_combined" => DataFrame),
    #  "classnames" -> list[str],
    #  "class_abundance -> dict[class_id, dict(CA:float, count: int)]
    # )
    results: dict[ConditionId, dict[str, Any]] = {}
    #: Pairwise comparisons of conditions
    # (condition1, condition2) => dict(
    #  "intersection_data",
    #  "metrics",
    #  "RLS_results", (RLS = Relocalization Score)
    #  "RLS_null",
    #  "nRLS_results",
    #  "nRLS_null",
    #  )
    comparison: dict[
        tuple[ConditionId, ConditionId], dict[str, pd.Series | pd.DataFrame]
    ] = {}
    #: Indicates which of the individual analysis steps
    #  have already been performed or not
    status: SessionStatusModel = SessionStatusModel()

    def reset_global_changes(self):
        self.comparison = {}
        self.status.comparison_global = False
        self.status.comparison_class = False

    def reset_static_statistics(self):
        self.reset_global_changes()
        self.results = {}
        self.status.proteome_prediction = False
        self.status.lipidome_prediction = False

    def reset_input_tp(self):
        self.tp_paths = []
        self.tp_tables = {}
        self.tp_pos = {}
        self.tp_data = {}

    def reset_input_fract(self):
        self.fract_paths = []
        self.fract_tables = {}
        self.fract_pos = {}
        self.fract_data = {}

    def reset_infract(self):
        self.fract_indata = {}
        self.fract_identifiers = {}

    def reset_intp(self):
        self.tp_indata = {}
        self.tp_identifiers = {}

    def reset_fract(self):
        self.fract_data = {"class": {}, "vis": {}}
        self.fract_std = {"class": {}, "vis": {}}
        self.fract_intermediate = {}
        self.fract_info = {}
        self.fract_conditions = []

    def reset_tp(self):
        self.tp_data = {}
        self.tp_intermediate = {}
        self.tp_info = pd.DataFrame()
        self.tp_conditions = []
        self.tp_icorr = {}

    def reset_fractionation(self):
        self.reset_fract()
        self.reset_marker()
        self.status.fractionation_data = False

    def reset_classification(self):
        self.reset_static_statistics()

        self.reset_global_changes()
        self.svm_marker = {}
        self.svm_test = {}
        self.svm_metrics = {}
        self.learning_xyz = {}

        self.status.training = False

    def reset_marker(self):
        self.marker_list = pd.DataFrame()
        self.fract_marker = {}
        self.fract_marker_vis = {}
        self.fract_test = {}
        self.fract_full = {}
        self.reset_classification()
        self.status.marker_matched = False

    def reset(self, other: SessionModel = None):
        """Reset to default values or copy from another session."""
        if other is None:
            other = SessionModel()

        for field_name, field_value in other:
            setattr(self, field_name, field_value)

    def to_numpy(self, filepath: Path | str):
        """Serialize using np.save."""
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        with open(filepath, "wb") as f:
            np.save(f, self.model_dump(), allow_pickle=True)

    @classmethod
    def from_numpy(cls, filepath: Path | str):
        """Deserialize using np.load."""
        filepath = Path(filepath)
        with open(filepath, "rb") as f:
            data = np.load(f, allow_pickle=True).item()
            return cls(**data)


def write_global_changes_reports(
    comparison: dict[
        tuple[ConditionId, ConditionId], dict[str, pd.Series | pd.DataFrame]
    ],
    outdir: Path | str,
) -> None:
    """Create Excel reports for the global changes."""
    for comb in comparison:
        fname = Path(
            outdir,
            f"CCMPS_comparison_{comb[0]}_{comb[1]}.xlsx",
        )
        selected_columns = [
            col
            for col in comparison[comb]["metrics"].columns
            if col.startswith("fRL_")
        ] + ["fRLS", "DS", "P(t)_RLS"]
        df_out = comparison[comb]["metrics"][selected_columns]
        df_out.columns = [
            col.replace("fRL_", "RL_Relocalization_")
            if col.startswith("fRL_")
            else "RLS_ReLocalizationScore"
            if col == "fRLS"
            else "DS_DistanceScore"
            if col == "DS"
            else "P-Value"
            if col == "P(t)_RLS"
            else col
            for col in df_out.columns
        ]
        df_out.to_excel(fname, index=True)


def write_class_changes_reports(
    model: SessionModel, outdir: Path | str
) -> None:
    """Create Excel reports for the class changes."""
    for condition in model.results:
        fname = Path(
            outdir,
            f"CCMPS_ClassComposition_{condition}.xlsx",
        )
        selected_columns = [
            col
            for col in model.results[condition]["metrics"].columns
            if col.startswith("nCPA")
        ] + ["TPA"]
        df_out = model.results[condition]["metrics"][selected_columns]
        df_out.columns = [
            col.replace(
                "nCPA_imp_",
                "nCPA_normalizedClasscentrigProteinAmount_",
            )
            if col.startswith("nCPA_")
            else "TPA_TotalProteinAmount"
            if col == "TPA"
            else col
            for col in df_out.columns
        ]
        df_out.to_excel(fname, index=True)

    for comb in model.comparison:
        fname = Path(
            outdir,
            f"CCMPS_ClassComparison_{comb[0]}_{comb[1]}.xlsx",
        )
        selected_columns = [
            col
            for col in model.comparison[comb]["metrics"].columns
            if col.startswith("nCFC_")
        ]
        df_out = model.comparison[comb]["metrics"][selected_columns]
        df_out.columns = [
            col.replace(
                "nCFC_",
                "nCFC_normalizedClasscentricFoldChange_",
            )
            if col.startswith("nCFC_")
            else col
            for col in df_out.columns
        ]
        df_out.to_excel(fname, index=True)


def write_comparison_reports(model: SessionModel, outdir: str | Path) -> None:
    for comb in model.comparison:
        fname = Path(
            outdir,
            f"CCMPS_comparison_{comb[0]}_{comb[1]}.tsv",
        )

        df_out = pd.DataFrame(
            index=model.comparison[comb]["intersection_data"].index
        )
        df_out = pd.merge(
            df_out,
            model.comparison[comb]["metrics"],
            left_index=True,
            right_index=True,
            how="left",
        )
        for colname in model.fract_info:
            df_out = pd.merge(
                df_out,
                model.fract_info[colname],
                left_index=True,
                right_index=True,
                how="left",
            )
        df_out.to_csv(
            fname,
            sep="\t",
            index=True,
            index_label="Identifier",
        )


def write_statistics_reports(model: SessionModel, outdir: str | Path) -> None:
    for condition in model.results:
        fname = Path(outdir, f"CCMPS_statistics_{condition}.tsv")
        df_out = pd.merge(
            model.fract_data["vis"][condition + "_median"],
            model.results[condition]["metrics"],
            left_index=True,
            right_index=True,
            how="outer",
        )
        for colname in model.fract_info:
            df_out = pd.merge(
                df_out,
                model.fract_info[colname],
                left_index=True,
                right_index=True,
                how="left",
            )
        df_out.to_csv(
            fname,
            sep="\t",
            index=True,
            index_label="Identifier",
        )
