"""Core classes and functions for the ccompass package."""

from __future__ import annotations

import copy
import logging
import tempfile
import uuid
import zipfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
import yaml
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_serializer,
)

from ._utils import get_ccmps_data_directory

logger = logging.getLogger(__name__)


#: the application settings file
config_filepath: Path = get_ccmps_data_directory() / "settings.yaml"
#: the repository URL
repository_url = "https://github.com/ICB-DCM/C-COMPASS/"
#: the ReadTheDocs URL
readthedocs_url = "https://c-compass.readthedocs.io/en/latest/"
#: name of the application
app_name = "C-COMPASS"


#: The value used in the sample tables as condition IDs to indicate the column
#  that contains the protein IDs.
IDENTIFIER = "[IDENTIFIER]"
#: The value used in the sample tables as condition IDs to indicate the columns
#  with values that should be carried forward to the final results.
#  (And potentially be used for matching markers.)
KEEP = "[KEEP]"
#: The value used in the sample table for not applicable values.
#  I.e. for "Replicate" and "Fraction" columns for KEEP and IDENTIFIER rows.
NA = "-"


class AppSettings(BaseModel):
    """Settings for the C-COMPASS application"""

    #: The directory that was last used to load/save a session
    last_session_dir: Path = Path.home()

    #: The maximum number of processes to use for parallel processing
    max_processes: int = 1

    #: Recently used files (sessions)
    recent_files: list[Path] = []

    #: Number of recent files to keep
    max_recent_files: int = 10

    def add_recent_file(self, filepath: Path | str):
        """Add a file to the list of recent files."""
        filepath = Path(filepath)
        if filepath in self.recent_files:
            self.recent_files.remove(filepath)
        self.recent_files.insert(0, filepath)
        self.recent_files = self.recent_files[: self.max_recent_files]

    @field_serializer("last_session_dir")
    def serialize_last_session_dir(self, value: Path) -> str:
        return str(value)

    @field_serializer("recent_files")
    def serialize_recent_files(self, value: list[Path]) -> list[str]:
        return list(map(str, value))

    @classmethod
    def load(cls, filepath: Path = None):
        """Load the settings from a file."""
        import yaml

        if filepath is None:
            filepath = config_filepath

        if not filepath.exists():
            return cls()

        logger.debug(f"Loading settings from {filepath}")

        with open(filepath) as f:
            data = yaml.safe_load(f) or {}
            return cls(**data)

    def save(self, filepath: Path = None):
        """Save the settings to a file."""
        import yaml

        if filepath is None:
            filepath = config_filepath

        filepath.parent.mkdir(parents=True, exist_ok=True)

        # create a backup of the old settings
        if filepath.exists():
            backup_path = filepath.with_suffix(".bak")
            backup_path.write_text(filepath.read_text())

        with open(filepath, "w") as f:
            yaml.safe_dump(self.model_dump(), f)


class NeuralNetworkParametersModel(BaseModel):
    """Hyperparameters for the neural network."""

    #: Perform upsampling?
    upsampling: bool = True
    #: Method for upsampling
    upsampling_method: Literal["none", "average", "noisedaverage"] = (
        "noisedaverage"
    )
    #: Noise level for upsampling (standard deviations)
    upsampling_noise: float = 2
    #: Perform SVM filtering?
    svm_filter: bool = False
    #: The number of different ratios for pairwise mixing,
    #  or "none" for no mixing.
    #  The ratios will be 1/N, 2/N, ..., (N-1)/N.
    mixed_part: int | str = 4
    #: The fraction of the mixed batch to use (∈ [0, 1])
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
    #: Optimizers to include in the hyperparameter search
    optimizers: list[Literal["adam", "rmsprop", "sgd"]] = [
        "adam",
        "rmsprop",
        "sgd",
    ]
    #: Number of epochs for the neural network training
    NN_epochs: int = 20
    #: The number of independent rounds for
    #  upsampling/mixing/training/prediction
    rounds: int = 3
    #: Repetitions for the neural network training to generate an ensemble
    subrounds: int = 10
    #: Percentile threshold for false-positive class probabilities
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


class XYZ_Model(BaseModel):
    """`learning_xyz` for a specific condition in `SessionModel`.

    W, Y: true labels
    w, y: predicted labels
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    condition_id: str = ""
    # TODO(performance): get rid of duplicate data in different formats
    #: List of unique classes for which there are marker measurements
    classes: list[str] = []
    #: class labels for the markers
    W_train_df: pd.Series = pd.Series()
    #: Class labels for all proteins (NaN for non-marker proteins)
    W_full_df: pd.Series = pd.Series()

    #: Combined classification results from different SVM rounds
    #  (w_full_prob_df)
    w_full_combined: pd.DataFrame = pd.DataFrame()
    #: Probabilities for the classifications in w_full_combined
    w_full_prob_combined: pd.DataFrame = pd.DataFrame()

    #: features (protein levels in the different fractions for one replicate,
    #  for proteins with known and unknown class labels)
    x_full_df: pd.DataFrame = pd.DataFrame()

    #: Features for the proteins with known class labels
    x_test_df: pd.DataFrame = pd.DataFrame()
    #: Features for the training data (marker profiles)
    x_train_df: pd.DataFrame = pd.DataFrame()

    #: One-hot encoded labels for marker profiles
    Z_train_df: pd.DataFrame = pd.DataFrame()
    #: Means of the z_full values across the different rounds
    z_full_mean_df: pd.DataFrame = pd.DataFrame()

    #: Results / intermediate data for the different training rounds
    #  round_id => TrainingRound_Model
    round_results: dict[str, TrainingRoundModel] = {}


class TrainingRoundModel(BaseModel):
    """Data for a single round of model training.

    A upsamling/mixing/training/prediction round for a single condition."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: features and SVM classification results
    # TODO: contains an unnecessary copy of x_full_df
    w_full_prob_df: pd.DataFrame = pd.DataFrame()

    #: Summary of the best neural network mode
    FNN_summary: str = ""
    #: Class labels for the upsampled training data
    #  (protein (index), 'class')
    W_train_up_df: pd.Series = pd.Series()
    #: Features for the upsampled full dataset
    #  (protein x fraction)
    x_full_up_df: pd.DataFrame = pd.DataFrame()
    #: Features of the upsampled training data
    x_train_up_df: pd.DataFrame = pd.DataFrame()
    #: Features for the training data (marker profiles) after maxing
    x_train_mixed_up_df: pd.DataFrame = pd.DataFrame()
    #: Class probabilities for mixed profiles (mixing ratios)
    #  (protein x compartment)
    Z_train_mixed_up_df: pd.DataFrame = pd.DataFrame()

    #: Data for the different rounds of neural network training after the
    #  hyperparameter search. Basis for ensemble prediction.
    subround_results: dict[str, TrainingSubRoundModel] = {}


class TrainingSubRoundModel(BaseModel):
    """Data for a single round of neural network model training and prediction."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: Neural network classification results for y_train
    #  (i.e. probabilities for the different classes for each protein)
    z_train_df: pd.DataFrame = pd.DataFrame()
    #: Neural network classification results for y_full
    #  (i.e. probabilities for the different classes for each protein)
    z_full_df: pd.DataFrame = pd.DataFrame()


class ResultsModel(BaseModel):
    """Results for a single condition (static statistics)."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: Table with:
    #  * `marker`: the compartment name if the row is a marker, NA otherwise
    #  * `SVM_winner`: the SVM-predicted compartment if SVM predictions
    #     agreed across replicates and rounds, NA otherwise
    #  * `SVM_prob`: the probability of the SVM-predicted compartment
    #  * `SVM_subwinner`: the compartment that was most often predicted
    #     by the different SVM rounds, NA in case of a tie
    #  * `CC_$class_$replicateId`: the multi-class predictions for the
    #       individual replicates
    #  * `CC_$class`: the mean across `CC_$class_$replicateId`
    #  * `NN_winner`: the class/compartment with the highest CC value
    #       according to the neural network predictions
    #  * `fCC_$class`: the false-positive filtered and renormalized
    #       class contribution
    #  * `fNN_winner`: the class/compartment with the highest fCC value
    metrics: pd.DataFrame = pd.DataFrame()
    class_abundance: dict[str, dict[str, float | int]] = {}
    # Unique list of class names
    classnames: list[str] = []
    # List of "subcondition" IDs ("{condition}_Rep.{replicate}")
    subcons: list[str] = []
    #: SVM results
    #  (the combined results of the different replicates and SVM rounds)
    #  * winner_combined: DataFrame
    #  * prob_combined: DataFrame
    SVM: dict[str, pd.DataFrame] = {}


class ComparisonModel(BaseModel):
    """Result of a comparison between two conditions."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    intersection_data: pd.DataFrame = pd.DataFrame()
    metrics: pd.DataFrame = pd.DataFrame()


class MarkerSet(BaseModel):
    """A single marker table with some ID and class annotations as provided by
    the user.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: the user-provided marker file
    df: pd.DataFrame = Field(
        default=pd.DataFrame(), serialization_alias="table"
    )
    #: column ID in `df` to match the fractionation data identifiers
    #  ("key column" in GUI)
    # "-" means unset
    identifier_col: str = NA
    #: column ID in `df` that contains the class names
    # "-" means unset
    class_col: str = NA

    @property
    def classes(self) -> list[str]:
        return self.df[self.class_col].unique().tolist()


class FractDataset(BaseModel):
    """A fractionation data table possibly containing multiple
    conditions, fractions and replicates.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    # TODO: consider moving GUI-specific data out of here
    #: List of rows for the fractionation table in the GUI, referring to the
    #  columns in `df`.
    #  [column ID (str), condition (str), replicate (int), fraction (int)]
    table: list[list[int | str]] = []
    #: The actual fractionation data table
    #  Rows are proteins, columns are samples.
    #  Additional columns may be present.
    df: pd.DataFrame = pd.DataFrame()

    @property
    def id_col(self) -> str:
        """The column ID that contains the protein identifiers."""
        for row in self.table:
            if row[1] == IDENTIFIER:
                return row[0]

    @id_col.setter
    def id_col(self, value: int | str):
        """Set the column ID that contains the protein identifiers."""
        # unset condition, fraction and replicate for previous identifier row
        prev_id_row = [
            i for i, row in enumerate(self.table) if row[1] == IDENTIFIER
        ]
        if prev_id_row:
            assert len(prev_id_row) == 1
            prev_id_row = prev_id_row[0]
            self.table[prev_id_row][1:4] = ["", "", ""]

        if isinstance(value, int):
            # set by index
            self.table[value][1:4] = [IDENTIFIER, NA, NA]
        else:
            # set by column ID
            for row in self.table:
                if row[0] == value:
                    row[1:4] = [IDENTIFIER, NA, NA]
                    break


class TotalProtDataset(BaseModel):
    """A total proteome data table possibly containing multiple conditions
    and replicates.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    # TODO: consider moving GUI-specific data out of here
    #: List of rows for the total proteome table in the GUI, referring to the
    #  columns in `df`.
    table: list[list[str]] = []
    #: The actual total proteome data table.
    #  Rows are proteins, columns are samples.
    #  Additional columns may be present.
    df: pd.DataFrame = pd.DataFrame()

    @property
    def id_col(self) -> str:
        """The column ID that contains the protein identifiers."""
        for row in self.table:
            if row[1] == IDENTIFIER:
                return row[0]

    @id_col.setter
    def id_col(self, value: int | str):
        """Set the column ID that contains the protein identifiers."""
        # unset condition for previous identifier row
        prev_id_row = [
            i for i, row in enumerate(self.table) if row[1] == IDENTIFIER
        ]
        if prev_id_row:
            assert len(prev_id_row) == 1
            prev_id_row = prev_id_row[0]
            self.table[prev_id_row][1] = ""

        if isinstance(value, int):
            # set by index
            self.table[value][1] = IDENTIFIER
        else:
            # set by column ID
            for row in self.table:
                if row[0] == value:
                    row[1] = IDENTIFIER
                    break


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

    ## User input fractionation data

    #: User-provided fractionation data tables
    fract_input: dict[Filepath, FractDataset] = {}

    #: Fractionation preprocessing parameters.
    #  global/classification/visualization
    #  "global"|"class"|"vis" => option => value
    fract_preparams: dict[str, dict[str, Any]] = fract_default()
    #: The column ID of the fractionation DataFrame that is
    #  to be used for matching the markers (`marker_list["name"])
    marker_fractkey: str = IDENTIFIER

    ## User input markers

    #: The user-provided marker sets
    marker_sets: dict[Filepath, MarkerSet] = {}
    #: Options for merging marker sets
    marker_params: dict[str, str] = {"how": "exclude", "what": "unite"}
    #: Mapping of compartment names to class names
    #  nan-values indicate that the compartment is not to be used
    marker_conv: dict[str, str | float] = {}

    ## Processed fractionation data

    #: Fractionation data for classification and visualization
    #  For classification, one DataFrame for each condition x replicate
    #  ("{condition}_Rep.{replicate}")
    #  For visualization, one DataFrame for each condition
    #  ("{condition}_median")
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
    #: Addition ("keep") columns from the fractionation data
    #  column ID => DataFrame
    fract_info: dict[str, pd.DataFrame] = {}
    #: Conditions in the fractionation data, including "[KEEP]"
    fract_conditions: list[str] = []

    #: Marker profiles for the different conditions / replicates
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_marker: dict[ConditionReplicate, pd.DataFrame] = {}
    #: Median of marker profiles across replicates for each condition
    #  for visualization
    #  "{condition}_median" => DataFrame
    fract_marker_vis: dict[str, pd.DataFrame] = {}
    #: Fractionation data for non-marker proteins
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_test: dict[ConditionReplicate, pd.DataFrame] = {}
    #: Combined marker and non-marker profiles for the different
    #  conditions x replicates
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_full: dict[ConditionReplicate, pd.DataFrame] = {}

    #: The consolidated marker list, after merging `marker_sets`
    #  according to `marker_params`, and accounting for renaming
    #  and filtering through `marker_conv`.
    #  "name" (gene name, index) => "class" (class name)
    marker_list: pd.DataFrame = pd.DataFrame()

    ## User input total proteome data

    #: User-provided total proteome tables
    tp_input: dict[Filepath, TotalProtDataset] = {}

    #: Total proteome preprocessing parameters
    tp_preparams: dict[str, Any] = {"minrep": 2, "imputation": "normal"}

    ## Processed total proteome data

    #: Total proteome data for the different conditions
    #  One DataFrame for each condition containing all replicates
    #  (column names are "{condition}_Rep.{replicate}")
    tp_data: dict[ConditionReplicate, pd.DataFrame] = {}
    #: ??
    tp_icorr: dict = {}
    #: The [KEEP] columns for the combined total proteome dataset
    tp_info: pd.DataFrame = pd.DataFrame()

    ## User input classification data
    #: Neural network hyperparameters
    NN_params: NeuralNetworkParametersModel = NeuralNetworkParametersModel()

    #: Neural network data
    # "{condition}_Rep.{replicate}" => dict(
    #  {w,W,x,X,y,Y,z,Z}_... => ...
    # )
    learning_xyz: dict[ConditionReplicate, XYZ_Model] = {}

    #: `stats_proteome` results for the different conditions
    results: dict[ConditionId, ResultsModel] = {}
    #: Pairwise comparisons of conditions
    # (condition1, condition2) => ComparisonModel
    comparison: dict[tuple[ConditionId, ConditionId], ComparisonModel] = {}

    @property
    def status(self) -> SessionStatusModel:
        """Return a status object that keeps track of the analysis steps."""
        return SessionStatusModel(
            fractionation_data=bool(self.fract_data["class"]),
            tp_data=bool(self.tp_data),
            lipidome_data=False,
            lipidome_total=False,
            marker_file=bool(self.marker_sets),
            marker_matched=bool(self.fract_full),
            training=bool(self.learning_xyz),
            proteome_prediction=bool(self.results),
            lipidome_prediction=False,
            comparison_global=bool(self.comparison),
            comparison_class=bool(
                all("TPA" in r.metrics for r in self.results.values())
            ),
        )

    @property
    def fract_paths(self) -> list[str]:
        """Return the paths of the fractionation data files."""
        return list(self.fract_input.keys())

    @property
    def tp_paths(self) -> list[str]:
        """Return the paths of the total proteome data files."""
        return list(self.tp_input.keys())

    def reset_class_centric_changes(self):
        """Reset class-centric analysis results."""
        results = self.results
        comparisons = self.comparison

        for condition, result in results.items():
            result.class_abundance = {}
            result.metrics.drop(["TPA", "CA_relevant"], axis=1, inplace=True)

            for classname in result.classnames:
                result.metrics.drop(
                    [
                        *[
                            f"nCC_{classname}_{subcon}"
                            for subcon in result.subcons
                        ],
                        "nCC_" + classname,
                        "CPA_" + classname,
                        "CPA_log_" + classname,
                        "CPA_imp_" + classname,
                        "nCPA_" + classname,
                        "nCPA_log_" + classname,
                        "nCPA_imp_" + classname,
                    ],
                    axis=1,
                    inplace=True,
                )

        for comb, comparison in comparisons.items():
            comparison.metrics.drop(
                ["nRLS", "P(t)_nRLS", "P(u)_nRLS"], axis=1, inplace=True
            )

            for classname in results[comb[0]].classnames:
                comparison.metrics.drop(
                    [
                        "nRL_" + classname,
                        "CFC_" + classname,
                        "nCFC_" + classname,
                    ],
                    axis=1,
                    inplace=True,
                )

    def reset_global_changes(self):
        self.comparison = {}

    def reset_static_statistics(self):
        self.reset_global_changes()
        self.results = {}

    def reset_input_tp(self):
        self.tp_input = {}
        self.tp_data = {}

    def reset_input_fract(self):
        self.fract_input = {}
        self.fract_data = {}

    def reset_fract(self):
        self.fract_data = {"class": {}, "vis": {}}
        self.fract_std = {"class": {}, "vis": {}}
        self.fract_info = {}
        self.fract_conditions = []

    def reset_tp(self):
        self.tp_data = {}
        self.tp_info = pd.DataFrame()
        self.tp_icorr = {}

    def reset_fractionation(self):
        self.reset_fract()
        self.reset_marker()

    def reset_classification(self):
        self.reset_static_statistics()
        self.reset_global_changes()
        self.learning_xyz = {}

    def reset_marker(self):
        self.marker_list = pd.DataFrame()
        self.fract_marker = {}
        self.fract_marker_vis = {}
        self.fract_test = {}
        self.fract_full = {}
        self.reset_classification()

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

    def to_zip(self, filepath: Path | str):
        """Serialize the model to a zip file with YAML, TSV, and numpy files."""
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)

            def dataframe_representer(dumper, data: pd.DataFrame):
                """Custom YAML representer for pandas DataFrames."""
                if data.empty:
                    return dumper.represent_scalar("!pandas.DataFrame", "")

                file_id = str(uuid.uuid4())
                file_path = temp_dir / f"{file_id}.tsv"
                data.to_csv(file_path, sep="\t", index=True)
                return dumper.represent_scalar(
                    "!pandas.DataFrame", file_path.name
                )

            def ndarray_representer(dumper, data):
                """Custom YAML representer for numpy arrays."""
                file_id = str(uuid.uuid4())
                file_path = temp_dir / f"{file_id}.npy"
                np.save(file_path, data, allow_pickle=False)
                return dumper.represent_scalar(
                    "!numpy.ndarray", file_path.name
                )

            def series_representer(dumper, data):
                """Custom YAML representer for pandas Series."""
                file_id = str(uuid.uuid4())
                file_path = temp_dir / f"{file_id}.tsv"
                data.to_csv(file_path, sep="\t", index=True)
                return dumper.represent_scalar(
                    "!pandas.Series", file_path.name
                )

            def float64_representer(dumper, data):
                return dumper.represent_float(float(data))

            def tuple_representer(dumper, data):
                return dumper.represent_sequence("!tuple", data)

            yaml.add_representer(
                np.float64, float64_representer, Dumper=yaml.SafeDumper
            )
            yaml.add_representer(
                pd.DataFrame, dataframe_representer, Dumper=yaml.SafeDumper
            )
            yaml.add_representer(
                np.ndarray, ndarray_representer, Dumper=yaml.SafeDumper
            )
            yaml.add_representer(
                pd.Series, series_representer, Dumper=yaml.SafeDumper
            )
            yaml.add_representer(
                tuple, tuple_representer, Dumper=yaml.SafeDumper
            )

            with open(temp_dir / "session.yaml", "w") as f:
                yaml.safe_dump(self.model_dump(), f)

            with zipfile.ZipFile(
                filepath, "w", compression=zipfile.ZIP_DEFLATED
            ) as zipf:
                for item in temp_dir.iterdir():
                    zipf.write(item, item.name)

    @classmethod
    def from_zip(cls, filepath: Path | str):
        """Deserialize the model from a zip file with YAML and TSV files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)

            with zipfile.ZipFile(filepath, "r") as zipf:
                zipf.extractall(temp_dir)

            def dataframe_constructor(loader, node):
                """Custom YAML constructor for pandas DataFrames."""
                if not (filename := loader.construct_scalar(node)):
                    return pd.DataFrame()

                file_path = temp_dir / filename
                try:
                    return pd.read_csv(file_path, sep="\t", index_col=0)
                except pd.errors.EmptyDataError:
                    return pd.DataFrame()

            def ndarray_constructor(loader, node):
                """Custom YAML constructor for numpy arrays."""
                file_path = temp_dir / loader.construct_scalar(node)
                return np.load(file_path, allow_pickle=False)

            def series_constructor(loader, node):
                """Custom YAML constructor for pandas Series."""
                file_path = temp_dir / loader.construct_scalar(node)
                df = pd.read_csv(
                    file_path,
                    sep="\t",
                    index_col=0,
                    header=0,
                    float_precision="round_trip",
                )
                assert df.shape[1] == 1
                return df.iloc[:, 0]

            def tuple_constructor(loader, node):
                return tuple(loader.construct_sequence(node))

            yaml.add_constructor(
                "!pandas.DataFrame",
                dataframe_constructor,
                Loader=yaml.SafeLoader,
            )
            yaml.add_constructor(
                "!numpy.ndarray", ndarray_constructor, Loader=yaml.SafeLoader
            )
            yaml.add_constructor(
                "!pandas.Series", series_constructor, Loader=yaml.SafeLoader
            )
            yaml.add_constructor(
                "!tuple", tuple_constructor, Loader=yaml.SafeLoader
            )

            with open(temp_dir / "session.yaml") as f:
                data = yaml.safe_load(f)
                return cls(**data)


def write_global_changes_reports(
    comparison: dict[tuple[ConditionId, ConditionId], ComparisonModel],
    outdir: Path | str,
) -> None:
    """Create Excel reports for the global changes."""
    Path(outdir).mkdir(parents=True, exist_ok=True)

    for comb in comparison:
        fname = Path(
            outdir,
            f"CCMPS_comparison_{comb[0]}_{comb[1]}.xlsx",
        )
        selected_columns = [
            col
            for col in comparison[comb].metrics.columns
            if col.startswith("fRL_")
        ] + ["fRLS", "DS", "P(t)_RLS"]
        df_out = comparison[comb].metrics[selected_columns]
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
    Path(outdir).mkdir(parents=True, exist_ok=True)

    for condition, result in model.results.items():
        fname = Path(
            outdir,
            f"CCMPS_ClassComposition_{condition}.xlsx",
        )
        selected_columns = [
            col for col in result.metrics.columns if col.startswith("nCPA")
        ] + ["TPA"]
        df_out = result.metrics[selected_columns]
        df_out.columns = [
            col.replace(
                "nCPA_imp_",
                "nCPA_normalizedClasscentricProteinAmount_",
            )
            if col.startswith("nCPA_")
            else "TPA_TotalProteinAmount"
            if col == "TPA"
            else col
            for col in df_out.columns
        ]
        df_out.to_excel(fname, index=True)

    for (cond1, cond2), comp in model.comparison.items():
        fname = Path(
            outdir,
            f"CCMPS_ClassComparison_{cond1}_{cond2}.xlsx",
        )
        selected_columns = [
            col for col in comp.metrics.columns if col.startswith("nCFC_")
        ]
        df_out = comp.metrics[selected_columns]
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
    Path(outdir).mkdir(parents=True, exist_ok=True)

    for (cond1, cond2), comp in model.comparison.items():
        fname = Path(
            outdir,
            f"CCMPS_comparison_{cond1}_{cond2}.tsv",
        )

        df_out = pd.DataFrame(index=comp.intersection_data.index)
        df_out = pd.merge(
            df_out,
            comp.metrics,
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
    Path(outdir).mkdir(parents=True, exist_ok=True)

    for condition, result in model.results.items():
        fname = Path(outdir, f"CCMPS_statistics_{condition}.tsv")
        df_out = pd.merge(
            model.fract_data["vis"][condition + "_median"],
            result.metrics,
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


def create_identity_conversion(
    marker_sets: Iterable[MarkerSet],
) -> dict[str, str]:
    """Create dummy conversion dictionary for classes."""
    return {
        classname: classname
        for marker_set in marker_sets
        for classname in marker_set.classes
    }


def create_markerlist(
    marker_sets: dict[str, MarkerSet],
    marker_conv: dict[str, str | float],
    what: Literal["unite", "intersect"],
    how: Literal["majority", "exclude"],
) -> pd.DataFrame:
    """Create a uniform marker list from multiple marker sets.

    Create a uniform marker list from multiple marker sets, accounting for
    any filtering or renaming.

    :param marker_sets: Marker sets to unify.
    :param marker_conv: Mapping of class names in the marker sets to the
        unified class names. Usually an identity mapping.
        NaN values indicate that the marker should be excluded.
    :param what: Whether to return the union or intersection of the markers.
    :param how: How to resolve conflicting class assignments.
        "majority": Assign the most frequent class.
        "exclude": Exclude the marker from the list.

    :return: A DataFrame with the unified marker list with marker names as
        index ('name') and a 'class' column.
    """
    if what not in ["unite", "intersect"]:
        raise ValueError(f"Invalid 'what' parameter: {what}")
    if how not in ["majority", "exclude"]:
        raise ValueError(f"Invalid 'how' parameter: {how}")

    # handle ID conversion, normalize column names
    combined = pd.DataFrame(columns=["name"])
    counter = 1
    for marker_set in marker_sets.values():
        id_col = marker_set.identifier_col
        class_col = marker_set.class_col
        cur_df = marker_set.df[[id_col, class_col]].copy()
        cur_df.rename(
            columns={
                id_col: "name",
                class_col: f"class{counter}",
            },
            inplace=True,
        )
        class_col = f"class{counter}"
        cur_df.replace(
            {class_col: marker_conv},
            inplace=True,
        )
        cur_df.replace(
            {class_col: {r"^\s*$": np.nan}},
            regex=True,
            inplace=True,
        )
        cur_df = cur_df[cur_df[class_col].notna()]

        combined = pd.merge(combined, cur_df, on="name", how="outer")
        counter += 1

    combined.set_index("name", inplace=True)

    # union or intersection of markers?
    if what == "unite":
        pass
    elif what == "intersect":
        combined.dropna(inplace=True)

    # resolve conflicting class assignments
    # majority: assign the most frequent class
    # exclude: exclude the marker from the list
    if how == "majority":
        combined = pd.DataFrame(combined.mode(axis=1, dropna=True)[0]).rename(
            columns={0: "class"}
        )
    elif how == "exclude":
        combined = combined.mode(axis=1, dropna=True).fillna(np.nan)
        if 1 in combined.columns:
            combined = pd.DataFrame(combined[combined[1].isnull()][0]).rename(
                columns={0: "class"}
            )
        else:
            combined.rename(columns={0: "class"}, inplace=True)

    return combined


def create_marker_profiles(
    fract_data,
    key: str,
    fract_info: dict[str, pd.DataFrame],
    marker_list: pd.DataFrame,
) -> tuple[
    dict[str, pd.DataFrame], dict[str, pd.DataFrame], dict[str, pd.DataFrame]
]:
    """Create marker profiles for classification and visualization."""
    # create marker profiles for classification
    profiles: dict[str, pd.DataFrame] = {}
    for condition in fract_data["class"]:
        profiles[condition] = copy.deepcopy(fract_data["class"][condition])

    # profiles with known class
    fract_marker: dict[str, pd.DataFrame] = {}
    # profiles with unknown class
    fract_test: dict[str, pd.DataFrame] = {}

    for condition in profiles:
        if key == IDENTIFIER:
            # match on index
            profile_full = pd.merge(
                profiles[condition],
                marker_list,
                left_index=True,
                right_index=True,
                how="left",
            )
        else:
            # match on some [KEEP] column from fract_info
            # add identifier column
            profiles[condition] = pd.merge(
                profiles[condition],
                fract_info[key].astype(str).map(str.upper),
                left_index=True,
                right_index=True,
            )
            profile_full = pd.merge(
                profiles[condition],
                marker_list,
                left_on=key,
                right_index=True,
                how="left",
            ).drop(key, axis=1)
        fract_marker[condition] = profile_full.dropna(subset=["class"])
        fract_test[condition] = profile_full[profile_full["class"].isna()]
        if fract_marker[condition].empty:
            raise ValueError(
                f"No markers found for condition {condition}. "
                "Please check the marker list, the selected keys, "
                "and the fractionation data."
            )
    # create marker profiles for visualization
    profiles_vis = {}
    for condition in fract_data["vis"]:
        profiles_vis[condition] = copy.deepcopy(fract_data["vis"][condition])

    fract_marker_vis: dict[str, pd.DataFrame] = {}
    for condition in profiles_vis:
        if key == IDENTIFIER:
            fract_marker_vis[condition] = pd.merge(
                profiles_vis[condition],
                marker_list,
                left_index=True,
                right_index=True,
            )
        else:
            profiles_vis[condition] = pd.merge(
                profiles_vis[condition],
                fract_info[key],
                left_index=True,
                right_index=True,
            )
            fract_marker_vis[condition] = (
                pd.merge(
                    profiles_vis[condition],
                    marker_list,
                    left_on=key,
                    right_index=True,
                    how="left",
                )
                .drop(key, axis=1)
                .dropna(subset=["class"])
            )

    return fract_marker, fract_marker_vis, fract_test


def create_fullprofiles(
    fract_marker: dict[str, pd.DataFrame], fract_test: dict[str, pd.DataFrame]
) -> dict[str, pd.DataFrame]:
    """Create full profiles for classification.

    Concatenate the marker profiles with the test profiles.
    """
    fract_full = {}
    for condition in fract_test:
        fract_full[condition] = pd.concat(
            [fract_test[condition], fract_marker[condition]]
        )
    return fract_full


def read_marker_table(
    filepath: Path | str,
) -> pd.DataFrame:
    """Read a marker table in tsv format."""
    df = pd.read_csv(filepath, sep="\t", header=0).apply(
        lambda x: x.astype(str).str.upper()
    )
    logger.debug(f"Read marker table from {filepath} ({df.shape})")
    return df


def read_fract_table(
    filepath: Path | str,
) -> pd.DataFrame:
    """Read a fractionation table in tsv format."""
    df = pd.read_csv(filepath, sep="\t", header=0)
    df = df.replace("NaN", np.nan)
    df = df.replace("Filtered", np.nan)
    logger.debug(f"Read fractionation data from {filepath} ({df.shape})")
    return df


def read_tp_table(
    filepath: Path | str,
) -> pd.DataFrame:
    """Read a total proteome table in tsv format."""
    df = pd.read_csv(filepath, sep="\t", header=0)
    df = df.replace("NaN", np.nan)
    df = df.replace("Filtered", np.nan)
    df = df.map(convert_to_float)
    rows_with_float = df.map(is_float).any(axis=1)
    df = df[rows_with_float]
    logger.debug(f"Read total proteome data from {filepath} ({df.shape})")
    return df


def is_float(element):
    try:
        float(element)
        return True
    except (ValueError, TypeError):
        return False


def convert_to_float(x):
    try:
        return float(x)
    except ValueError:
        return x
