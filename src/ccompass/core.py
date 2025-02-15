"""Core classes and functions for the ccompass package."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, field_serializer, field_validator

from . import config_filepath

logger = logging.getLogger(__name__)


class AppSettings(BaseModel):
    """Settings for the C-COMPASS application"""

    #: The directory that was last used to load/save a session
    last_session_dir: Path = Path.home()

    #: The maximum number of processes to use for parallel processing
    max_processes: int = 1

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
    #: Repetitions for the neural network training to generate an ensemble
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


class XYZ_Model(BaseModel):
    """`learning_xyz` for a specific condition in `SessionModel`.

    W, Y: true labels
    w, y: predicted labels
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    # TODO(performance): get rid of duplicate data in different formats
    #: List of unique classes for which there are marker measurements
    classes: list[str] = []
    #: class labels for the markers
    W_train_df: pd.Series = pd.Series()
    #: class labels for the markers as list
    W_train: list = []
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
    #: x_full_df, but as numpy array
    x_full: np.ndarray = np.array([])

    #: Features for the proteins with known class labels
    x_test_df: pd.DataFrame = pd.DataFrame()
    #: x_test_df as numpy array
    x_test: np.ndarray = np.array([])
    #: Features for the training data (marker profiles)
    x_train_df: pd.DataFrame = pd.DataFrame()
    #: x_train_df as numpy array
    x_train: np.ndarray = np.array([])

    #: One-hot encoded labels for marker profiles
    Z_train_df: pd.DataFrame = pd.DataFrame()
    #: Means of the z_full values across the different rounds
    z_full_mean_df: pd.DataFrame = pd.DataFrame()

    #: Results / intermediate data for the different training rounds
    #  round_id => TrainingRound_Model
    round_results: dict[str, TrainingRound_Model] = {}


class TrainingRound_Model(BaseModel):
    """Data for a single round of model training.

    A upsamling/mixing/training/prediction round for a single condition."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: features and SVM classification results
    w_full_prob_df: pd.DataFrame = pd.DataFrame()

    #: Summary of the best neural network mode
    FNN_summary: str = ""
    #: Class labels for the upsampled training data
    W_train_up_df: pd.Series = pd.Series()
    #: W_train_up_df as list
    W_train_up: list = []
    #: Features for the upsampled full dataset
    x_full_up_df: pd.DataFrame = pd.DataFrame()
    #: x_full_up_df as numpy array
    x_full_up: np.ndarray = np.array([])
    #: Features of the upsampled training data
    x_train_up_df: pd.DataFrame = pd.DataFrame()
    #: x_train_up_df as numpy array
    x_train_up: np.ndarray = np.array([])
    #: same as x_full_up
    V_full_up: np.ndarray = np.array([])
    #: Features for the training data (marker profiles) after maxing
    x_train_mixed_up_df: pd.DataFrame = pd.DataFrame()
    #: x_train_mixed_up_df as numpy
    x_train_mixed_up: np.ndarray = np.array([])
    #: Class probabilities for mixed profiles (mixing ratios)
    Z_train_mixed_up_df: pd.DataFrame = pd.DataFrame()
    #: Z_train_mixed_up_df as numpy array
    Z_train_mixed_up: np.ndarray = np.array([])

    #: SVM-predicted class labels for x_full for each round
    # TODO: never read
    w_full: list = []
    #: Probabilities for the SVM-predicted class labels in w_full
    # TODO: never read
    w_full_prob: np.ndarray = np.array([])
    #: SVM-predicted class labels for x_train for each round
    # TODO: never read
    w_train: list = []

    #: Data for the different rounds of neural network training after the
    #  hyperparameter search. Basis for ensemble prediction.
    subround_results: dict[str, TrainingSubRound_Model] = {}

    @field_validator("w_full_prob", mode="before")
    def convert_lists_to_arrays(cls, v):
        # for backward compatibility, convert list to ndarray
        return {k: np.array(vv) for k, vv in v.items()}


class TrainingSubRound_Model(BaseModel):
    """Data for a single round of neural network model training and prediction."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    #: Class labels for the mixed profiles
    y_full_df: pd.DataFrame = pd.DataFrame()
    #: same as x_full
    y_full: np.ndarray = np.array([])
    #: same as x_full_up
    y_full_up: np.ndarray = np.array([])
    #: same as x_train_df
    y_train_df: pd.DataFrame = pd.DataFrame()
    #: same as x_train
    y_train: np.ndarray = np.array([])
    #: same as x_train_up
    y_train_up: np.ndarray = np.array([])
    #: same as x_train_mixed_up
    y_train_mixed_up: np.ndarray = np.array([])
    #: same as x_test
    y_test: np.ndarray = np.array([])

    #: Neural network classification results for y_train
    #  (i.e. probabilities for the different classes for each protein)
    z_train_df: pd.DataFrame = pd.DataFrame()
    #: Neural network classification results for y_full
    #  (i.e. probabilities for the different classes for each protein)
    z_full_df: pd.DataFrame = pd.DataFrame()
    #: z_full_df, but as numpy array
    z_full: np.ndarray = np.array([])


class ResultsModel(BaseModel):
    """Results for a single condition."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    metrics: pd.DataFrame = pd.DataFrame()
    class_abundance: dict[str, dict[str, dict[str, Any]]] = {}
    classnames: list[str] = []
    #: SVM results
    # * winner_combined: DataFrame
    # * prob_combined: DataFrame
    SVM: dict[str, pd.DataFrame] = {}


class ComparisonModel(BaseModel):
    """Result of a comparison between two conditions."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    intersection_data: pd.DataFrame = pd.DataFrame()
    metrics: pd.DataFrame = pd.DataFrame()

    # global comparison results
    #: Relocalization scores
    RLS_results: pd.Series = pd.Series()
    RLS_null: pd.Series = pd.Series()

    # class-centric comparison results
    nRLS_results: pd.Series = pd.Series()
    nRLS_null: pd.Series = pd.Series()


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

    # TODO: consider grouping data by condition

    ## User input fractionation data

    #: Fractionation column assignments
    #  filepath => [column ID, condition, replicate, fraction]
    fract_tables: dict[Filepath, list[list[int | str]]] = {}
    #: ??
    fract_pos: dict[Filepath, list[int]] = {}
    #: Fractionation input files: filepath => DataFrame
    fract_indata: dict[Filepath, pd.DataFrame] = {}
    #: Identifier column of each fractionation data table:
    #   filepath => column id
    fract_identifiers: dict[Filepath, str] = {}
    #: Fractionation preprocessing parameters.
    #  global/classification/visualization
    #  "global"|"class"|"vis" => option => value
    fract_preparams: dict[str, dict[str, Any]] = fract_default()
    #: The column ID of the fractionation DataFrame that is
    #  to be used for matching the markers (`marker_list["name"])
    marker_fractkey: str = "[IDENTIFIER]"

    ## User input markers

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

    ## Processed fractionation data

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
    #: Addition ("keep") columns from the fractionation data
    #  column ID => DataFrame
    fract_info: dict[str, pd.DataFrame] = {}
    #: Conditions in the fractionation data, including "[KEEP]"
    fract_conditions: list[str] = []

    #: Fractionation data for the different conditions x replicates
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_full: dict[ConditionReplicate, pd.DataFrame] = {}
    #: Marker abundance in the different fractions
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_marker: dict[ConditionReplicate, pd.DataFrame] = {}
    #: Marker abundance in the different fractions for visualization
    #  "{condition}_median" => DataFrame
    fract_marker_vis: dict[str, pd.DataFrame] = {}

    #: ??
    #  "{condition}_Rep.{replicate}" => DataFrame
    fract_test: dict[ConditionReplicate, pd.DataFrame] = {}

    #: The consolidated marker list, after merging `marker_sets`
    #  according to `marker_params`, and accounting for renaming
    #  and filtering through `marker_conv`.
    #  "name" (gene name, index) => "class" (class name)
    marker_list: pd.DataFrame = pd.DataFrame()

    ## User input total proteome data

    #: Filepaths for total proteome data
    tp_paths: list[Filepath] = []
    #: Total proteome column assignments
    #  filepath => [column ID, condition]
    tp_tables: dict[Filepath, list[list[str]]] = {}
    #: ??
    tp_pos: dict[Filepath, list[int]] = {}
    #: Total proteome input files: filepath => DataFrame
    tp_indata: dict[Filepath, pd.DataFrame] = {}
    #: Identifier column for the total proteome: filepath => column id
    tp_identifiers: dict[Filepath, str] = {}
    #: Total proteome preprocessing parameters
    tp_preparams: dict[str, Any] = {"minrep": 2, "imputation": "normal"}

    ## Processed total proteome data

    #: Total proteome data for the different conditions
    #  One DataFrame for each condition containing all replicates
    #  (column names are "{condition}_Rep.{replicate}")
    tp_data: dict[ConditionReplicate, pd.DataFrame] = {}
    #: ??
    tp_icorr: dict = {}
    #: ??
    tp_conditions: list = []
    #: ??
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
        self.fract_tables = {}
        self.fract_pos = {}
        self.fract_data = {}
        self.fract_indata = {}
        self.fract_identifiers = {}

    def reset_intp(self):
        self.tp_indata = {}
        self.tp_identifiers = {}

    def reset_fract(self):
        self.fract_data = {"class": {}, "vis": {}}
        self.fract_std = {"class": {}, "vis": {}}
        self.fract_info = {}
        self.fract_conditions = []

    def reset_tp(self):
        self.tp_data = {}
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
