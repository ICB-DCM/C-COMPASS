"""Multiple organelle prediction."""

import copy
import logging
from datetime import datetime

import keras_tuner as kt
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn import svm
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
)
from tensorflow import keras
from tensorflow.keras import ops
from tensorflow.keras.backend import epsilon

from ._utils import get_ccmps_data_directory
from .core import NeuralNetworkParametersModel

logger = logging.getLogger(__package__)


#: mapping of optimizer names to optimizer classes
optimizer_classes = {
    "adam": tf.keras.optimizers.Adam,
    "rmsprop": tf.keras.optimizers.RMSprop,
    "sgd": tf.keras.optimizers.SGD,
}


def _create_classifier_hypermodel(
    NN_params: NeuralNetworkParametersModel,
) -> type[kt.HyperModel]:
    """Create a hypermodel for the classifier."""

    class FNN_classifier(kt.HyperModel):
        def __init__(self, fixed_hp=None, set_shapes=None):
            super().__init__()
            self.fixed_hp = fixed_hp
            self.set_shapes = set_shapes
            self.chosen_hp = {}

        def build(self, hp):
            model = keras.Sequential()
            # Input layer, size is the number of fractions
            model.add(
                tf.keras.Input(
                    (self.set_shapes[0],),
                )
            )

            # fixed or tunable hyperparameters
            if self.fixed_hp:
                optimizer_choice = self.fixed_hp["optimizer"]
                learning_rate = self.fixed_hp["learning_rate"]
                units = self.fixed_hp["units"]
            else:
                optimizer_choice = hp.Choice("optimizer", NN_params.optimizers)
                learning_rate = hp.Float(
                    "learning_rate",
                    min_value=1e-4,
                    max_value=1e-1,
                    sampling="log",
                )
                if NN_params.NN_optimization == "short":
                    units = hp.Int(
                        "units",
                        min_value=int(
                            min(self.set_shapes)
                            + 0.4
                            * (max(self.set_shapes) - min(self.set_shapes))
                        ),
                        max_value=int(
                            min(self.set_shapes)
                            + 0.6
                            * (max(self.set_shapes) - min(self.set_shapes))
                        ),
                        step=2,
                    )
                elif NN_params.NN_optimization == "long":
                    units = hp.Int(
                        "units",
                        min_value=min(self.set_shapes),
                        max_value=max(self.set_shapes),
                        step=2,
                    )
                else:
                    raise ValueError(
                        f"Unknown optimization: {NN_params.NN_optimization}"
                    )

            # dense layer 1 with tunable size
            if NN_params.NN_activation == "relu":
                model.add(keras.layers.Dense(units, activation="relu"))
            elif NN_params.NN_activation == "leakyrelu":
                hp_alpha = hp.Float(
                    "alpha", min_value=0.05, max_value=0.3, step=0.05
                )
                model.add(keras.layers.Dense(units))
                model.add(keras.layers.LeakyReLU(hp_alpha))

            # dense layer 2 with size according to the number of compartments
            model.add(
                keras.layers.Dense(
                    self.set_shapes[1],
                    activation=NN_params.class_activation,
                )
            )
            model.add(keras.layers.ReLU())

            # normalization layer
            model.add(keras.layers.Lambda(sum1_normalization))

            optimizer = optimizer_classes[optimizer_choice](
                learning_rate=learning_rate
            )
            model.compile(
                loss=NN_params.class_loss,
                optimizer=optimizer,
                metrics=[
                    tf.keras.metrics.MeanSquaredError(),
                    tf.keras.metrics.MeanAbsoluteError(),
                ],
            )

            if not self.fixed_hp:
                self.chosen_hp = {
                    "optimizer": optimizer_choice,
                    "learning_rate": learning_rate,
                    "units": units,
                }

            return model

        def get_chosen_hyperparameters(self):
            return self.chosen_hp

    return FNN_classifier


def upsample_condition(
    NN_params: NeuralNetworkParametersModel,
    stds: pd.DataFrame,
    fract_full: pd.DataFrame,
    fract_marker: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Perform upsampling for the given condition."""
    fract_full_up = fract_full
    fract_marker_up = fract_marker

    class_sizes = fract_marker["class"].value_counts()
    class_maxsize = class_sizes.max()
    k = 1
    for classname, data_class in fract_marker.groupby("class"):
        if not (class_difference := class_maxsize - class_sizes[classname]):
            continue

        data_class = data_class.drop(columns=["class"])
        class_up = pd.DataFrame(columns=data_class.columns)
        # TODO only compute where necessary
        class_std = data_class.std(axis=0).to_frame().transpose()
        class_std_flat = class_std.values.flatten()

        for i in range(class_difference):
            if NN_params.upsampling_method == "noised":
                sample = data_class.sample(n=1)
                ID_rnd = sample.index[0]
                name_up = f"up_{k}_{ID_rnd}"
                k += 1

                profile_rnd_flat = sample.values.flatten()
                std_rnd = stds.loc[[ID_rnd]]
                std_rnd = std_rnd[~std_rnd.index.duplicated(keep="first")]
                std_rnd_flat = std_rnd.values.flatten()
                std_rnd_flat = np.tile(
                    std_rnd_flat,
                    int(profile_rnd_flat.size / std_rnd_flat.size),
                )

                nv = np.random.normal(
                    profile_rnd_flat,
                    NN_params.upsampling_noise * std_rnd_flat,
                    size=sample.shape,
                )
                nv = np.where(nv > 1, 1, np.where(nv < 0, 0, nv))
                profile_up = pd.DataFrame(nv, columns=sample.columns)

            elif NN_params.upsampling_method == "average":
                sample = data_class.sample(n=3, replace=True)
                name_up = f"up_{k}_{'_'.join(sample.index)}"
                k += 1
                profile_up = sample.median(axis=0).to_frame().transpose()

            elif NN_params.upsampling_method == "noisedaverage":
                sample = data_class.sample(n=3, replace=True)
                name_up = f"up_{k}_{'_'.join(sample.index)}"
                k += 1

                profile_av = sample.median(axis=0).to_frame().transpose()
                profile_av_flat = profile_av.values.flatten()
                nv = np.random.normal(
                    profile_av_flat,
                    NN_params.upsampling_noise * class_std_flat,
                    size=profile_av.shape,
                )
                nv = np.where(nv > 1, 1, np.where(nv < 0, 0, nv))
                profile_up = pd.DataFrame(nv, columns=profile_av.columns)
            else:
                raise ValueError(
                    f"Unknown upsampling method: {NN_params.upsampling_method}"
                )

            profile_up.index = [name_up]
            profile_up["class"] = [classname]
            assert len(profile_up) == 1
            # TODO(performance) use concat only once
            class_up = (
                pd.concat([class_up, profile_up], axis=0, ignore_index=False)
                if not class_up.empty
                else profile_up
            )

        fract_marker_up = pd.concat(
            [fract_marker_up, class_up],
            axis=0,
            ignore_index=False,
        )
        fract_full_up = pd.concat(
            [fract_full_up, class_up],
            axis=0,
            ignore_index=False,
        )
    # TODO: seems unnecessary?!
    fract_marker_up = fract_marker_up.sample(frac=1)
    fract_full_up = fract_full_up.sample(frac=1)

    assert len(fract_marker_up["class"].value_counts().unique()) == 1
    assert len(fract_full_up["class"].value_counts().unique()) == 1

    return fract_marker_up, fract_full_up


def sum1_normalization(x):
    """Normalize the input to sum to 1."""
    return x / (ops.sum(x, axis=1, keepdims=True) + epsilon())


def init_learning_xyz(conditions: list[str]) -> dict[str, dict[str, dict]]:
    learning_xyz = {}
    for condition in conditions:
        learning_xyz[condition] = {}
        learning_xyz[condition]["W_full_up_df"] = {}
        learning_xyz[condition]["W_full_up"] = {}
        learning_xyz[condition]["W_train_up_df"] = {}
        learning_xyz[condition]["W_train_up"] = {}
        learning_xyz[condition]["w_full"] = {}
        learning_xyz[condition]["w_full_prob"] = {}
        learning_xyz[condition]["w_full_prob_df"] = {}
        learning_xyz[condition]["w_train"] = {}
        learning_xyz[condition]["w_train_prob"] = {}
        learning_xyz[condition]["w_test"] = {}
        learning_xyz[condition]["w_test_prob"] = {}
        learning_xyz[condition]["x_full_up_df"] = {}
        learning_xyz[condition]["x_full_up"] = {}
        learning_xyz[condition]["x_train_up_df"] = {}
        learning_xyz[condition]["x_train_up"] = {}
        learning_xyz[condition]["Z_train_df"] = {}
        learning_xyz[condition]["Z_train"] = {}
        learning_xyz[condition]["Z_train_up"] = {}
        learning_xyz[condition]["V_full_up"] = {}
        learning_xyz[condition]["x_train_mixed_up_df"] = {}
        learning_xyz[condition]["x_train_mixed_up"] = {}
        learning_xyz[condition]["Z_train_mixed_up_df"] = {}
        learning_xyz[condition]["Z_train_mixed_up"] = {}
        learning_xyz[condition]["AE_summary"] = {}
        learning_xyz[condition]["AE_history"] = {}
        learning_xyz[condition]["y_full_df"] = {}
        learning_xyz[condition]["y_full"] = {}
        learning_xyz[condition]["y_full_up"] = {}
        learning_xyz[condition]["y_train_df"] = {}
        learning_xyz[condition]["y_train"] = {}
        learning_xyz[condition]["y_train_up"] = {}
        learning_xyz[condition]["y_train_mixed_up"] = {}
        learning_xyz[condition]["y_test"] = {}
        learning_xyz[condition]["FNN_summary"] = {}
        learning_xyz[condition]["FNN_history"] = {}
        learning_xyz[condition]["z_full_df"] = {}
        learning_xyz[condition]["z_full"] = {}
        learning_xyz[condition]["z_train_df"] = {}
        learning_xyz[condition]["z_train"] = {}

    return learning_xyz


def MOP_exec(
    fract_full: dict[str, pd.DataFrame],
    fract_marker_old: dict[str, pd.DataFrame],
    fract_test: dict[str, pd.DataFrame],
    stds: dict[str, pd.DataFrame],
    NN_params: NeuralNetworkParametersModel,
):
    """Perform multi-organelle prediction.

    :param fract_full: dictionary of full profiles
    """
    conditions = list(fract_full.keys())
    learning_xyz = init_learning_xyz(conditions)

    for i_round in range(1, NN_params.rounds + 1):
        logger.info(f"Executing round {i_round}...")

        fract_full_up = {}
        fract_marker_up = {}

        fract_marker = copy.deepcopy(fract_marker_old)

        if NN_params.upsampling:
            # upsample fractionation data for each condition x replicate
            for condition in conditions:
                logger.info(f"Upsampling condition {condition}")
                fract_marker_up[condition], fract_full_up[condition] = (
                    upsample_condition(
                        NN_params,
                        stds.get(condition),
                        fract_full[condition],
                        fract_marker[condition],
                    )
                )
            logger.info("upsampling done!")
        else:
            fract_marker_up = copy.deepcopy(fract_marker)
            fract_full_up = copy.deepcopy(fract_full)

        logger.info("creating data...")
        for condition in conditions:
            learning_xyz, classes = create_learninglist(
                learning_xyz,
                fract_full,
                fract_full_up,
                fract_marker,
                fract_marker_up,
                fract_test,
                condition,
                i_round,
            )
        logger.info("data complete!")

        svm_metrics = {}
        svm_marker = {}
        svm_test = {}
        for condition in conditions:
            clf = svm.SVC(kernel="rbf", probability=True)

            svm_metrics, svm_marker, svm_test = single_prediction(
                learning_xyz[condition],
                clf,
                svm_metrics,
                fract_marker,
                svm_marker,
                fract_test,
                svm_test,
                condition,
                i_round,
            )

        if NN_params.svm_filter:
            logger.info("Applying SVM filter...")
            # Remove the markers that are not predicted correctly by the SVM
            #  and upsample the rest
            fract_full_up = {}
            fract_marker_up = {}
            fract_marker_filtered = {}
            for condition in conditions:
                rows_to_keep = (
                    svm_marker[condition]["class"]
                    == svm_marker[condition]["svm_prediction"]
                )
                fract_marker_filtered[condition] = fract_marker[condition][
                    rows_to_keep
                ]
                logger.info(f"Upsampling condition {condition}")
                fract_marker_up[condition], fract_full_up[condition] = (
                    upsample_condition(
                        NN_params,
                        stds.get(condition),
                        fract_full[condition],
                        fract_marker_filtered[condition],
                    )
                )
            logger.info("SVM filtering done.")

        fract_unmixed_up = {}
        for condition in conditions:
            unmixed_dummies = pd.get_dummies(
                fract_marker_up[condition]["class"]
            )
            fract_unmixed_up[condition] = pd.concat(
                [
                    fract_marker_up[condition].drop("class", axis=1),
                    unmixed_dummies,
                ],
                axis=1,
            )

        if NN_params.mixed_part == "none":
            fract_mixed_up = copy.deepcopy(fract_unmixed_up)
        else:
            logger.info("Mixing profiles...")
            fract_mixed_up = {}
            mix_steps = [
                i / NN_params.mixed_part
                for i in range(1, NN_params.mixed_part)
            ]
            for condition in conditions:
                fract_mixed_up = mix_profiles(
                    mix_steps,
                    NN_params,
                    fract_marker_up,
                    fract_unmixed_up,
                    fract_mixed_up,
                    condition,
                )

        round_id = f"ROUND_{i_round}"
        for condition in conditions:
            learning_xyz[condition]["x_train_mixed_up_df"][round_id] = (
                fract_mixed_up[condition].drop(columns=classes)
            )
            learning_xyz[condition]["x_train_mixed_up"][round_id] = (
                learning_xyz[condition]["x_train_mixed_up_df"][
                    round_id
                ].to_numpy(dtype=float)
            )
            learning_xyz[condition]["Z_train_mixed_up_df"][round_id] = (
                fract_mixed_up[condition][classes]
            )
            learning_xyz[condition]["Z_train_mixed_up"][round_id] = (
                learning_xyz[condition]["Z_train_mixed_up_df"][
                    round_id
                ].to_numpy(dtype=float)
            )
        logger.info("mixing done!")

        for condition in conditions:
            x_full = learning_xyz[condition]["x_full"]
            x_full_up = learning_xyz[condition]["x_full_up"][round_id]
            x_train = learning_xyz[condition]["x_train"]
            x_train_up = learning_xyz[condition]["x_train_up"][round_id]
            x_train_mixed_up = learning_xyz[condition]["x_train_mixed_up"][
                round_id
            ]
            x_test = learning_xyz[condition]["x_test"]

            V_full_up = learning_xyz[condition]["x_full_up"][round_id]
            learning_xyz[condition]["V_full_up"][round_id] = V_full_up

            if NN_params.AE == "none":
                y_full = x_full
                y_full_up = x_full_up
                y_train = x_train
                y_train_up = x_train_up
                y_train_mixed_up = x_train_mixed_up
                y_test = x_test

                learning_xyz = add_Y(
                    learning_xyz,
                    y_full,
                    y_full_up,
                    y_train,
                    y_train_up,
                    y_train_mixed_up,
                    y_test,
                    condition,
                    i_round,
                    0,
                )
                for SR in range(1, NN_params.subrounds + 1):
                    learning_xyz = add_Y(
                        learning_xyz,
                        y_full,
                        y_full_up,
                        y_train,
                        y_train_up,
                        y_train_mixed_up,
                        y_test,
                        condition,
                        i_round,
                        SR,
                    )
            else:
                # TODO ADD AUTOENCODER HERE
                raise NotImplementedError("Autoencoder not implemented yet.")

        FNN_classifier = _create_classifier_hypermodel(NN_params)
        for condition in conditions:
            learning_xyz = multi_predictions(
                FNN_classifier, learning_xyz, NN_params, condition, i_round
            )
    # TODO: no need to return / store fract_mixed_up, fract_unmixed_up, ...
    #  they are not used anymore
    return (
        learning_xyz,
        fract_full_up,
        fract_marker_up,
        fract_mixed_up,
        fract_unmixed_up,
        svm_marker,
        svm_test,
        svm_metrics,
    )


def multi_predictions(
    FNN_classifier: type[kt.HyperModel],
    learning_xyz,
    NN_params: NeuralNetworkParametersModel,
    condition: str,
    roundn: int,
):
    """Perform multi organelle predictions."""
    logger.info(f"Training classifier for condition {condition}...")
    y_full = learning_xyz[condition]["y_full"][f"ROUND_{roundn}_0"]
    y_train = learning_xyz[condition]["y_train"][f"ROUND_{roundn}_0"]
    y_train_mixed_up = learning_xyz[condition]["y_train_mixed_up"][
        f"ROUND_{roundn}_0"
    ]

    Z_train_mixed_up = learning_xyz[condition]["Z_train_mixed_up"][
        f"ROUND_{roundn}"
    ]
    set_shapes = [np.shape(y_train_mixed_up)[1], np.shape(Z_train_mixed_up)[1]]

    # Tune the hyperparameters
    classifier_directory = get_ccmps_data_directory()
    classifier_directory.mkdir(exist_ok=True, parents=True)

    now = datetime.now()
    time = now.strftime("%Y%m%d%H%M%S")
    tuner = kt.Hyperband(
        hypermodel=FNN_classifier(set_shapes=set_shapes),
        hyperparameters=kt.HyperParameters(),
        objective="val_mean_squared_error",
        max_epochs=NN_params.NN_epochs,
        factor=3,
        directory=str(classifier_directory),
        project_name=f"{time}_Classifier_{condition}_{roundn}",
    )

    stop_early = tf.keras.callbacks.EarlyStopping(
        monitor="val_loss", patience=5
    )

    tuner.search(
        y_train_mixed_up,
        Z_train_mixed_up,
        epochs=NN_params.NN_epochs,
        validation_split=0.2,
        callbacks=[stop_early],
    )
    logger.info("Hyperparameter tuning done!")
    best_model = tuner.get_best_models(num_models=1)[0]
    best_hp = tuner.get_best_hyperparameters(num_trials=1)[0].values

    stringlist = []
    best_model.summary(print_fn=lambda x: stringlist.append(x))
    FNN_summary = "\n".join(stringlist)

    learning_xyz[condition]["FNN_summary"]["ROUND_roundn"] = FNN_summary

    z_full = best_model.predict(y_full)
    z_train = best_model.predict(y_train)

    learning_xyz = add_Z(
        learning_xyz,
        z_full,
        z_train,
        condition,
        roundn,
        0,
    )

    for subround in range(1, NN_params.subrounds + 1):
        logger.info(
            f"Training classifier for condition {condition} "
            f"{subround}/{NN_params.subrounds}..."
        )
        # print(learning_xyz[condition]['y_full'])
        subround_id = f"ROUND_{roundn}_{subround}"
        y_full = learning_xyz[condition]["y_full"][subround_id]
        y_train = learning_xyz[condition]["y_train"][subround_id]
        y_train_mixed_up = learning_xyz[condition]["y_train_mixed_up"][
            subround_id
        ]

        fixed_model = FNN_classifier(
            fixed_hp=best_hp, set_shapes=set_shapes
        ).build(None)
        fixed_model.fit(
            y_train_mixed_up,
            Z_train_mixed_up,
            epochs=NN_params.NN_epochs,
            validation_split=0.2,
            callbacks=[stop_early],
        )

        z_full = fixed_model.predict(y_full)
        z_train = fixed_model.predict(y_train)

        learning_xyz = add_Z(
            learning_xyz,
            z_full,
            z_train,
            condition,
            roundn,
            subround,
        )

    return learning_xyz


def add_Z(
    learning_xyz,
    z_full,
    z_train,
    condition,
    roundn: int,
    subroundn: int,
):
    subround_id = f"ROUND_{roundn}_{subroundn}"
    learning_xyz[condition]["z_full_df"][subround_id] = pd.DataFrame(
        z_full,
        index=learning_xyz[condition]["y_full_df"][subround_id].index,
        columns=learning_xyz[condition]["classes"],
    )
    learning_xyz[condition]["z_full"][subround_id] = z_full

    learning_xyz[condition]["z_train_df"][subround_id] = pd.DataFrame(
        z_train,
        index=learning_xyz[condition]["y_train_df"][subround_id].index,
        columns=learning_xyz[condition]["classes"],
    )
    learning_xyz[condition]["z_train"][subround_id] = z_train

    return learning_xyz


def add_Y(
    learning_xyz,
    y_full,
    y_full_up,
    y_train,
    y_train_up,
    y_train_mixed_up,
    y_test,
    condition,
    roundn: int,
    subroundn: int,
):
    subround_id = f"ROUND_{roundn}_{subroundn}"
    learning_xyz[condition]["y_full_df"][subround_id] = pd.DataFrame(
        y_full, index=learning_xyz[condition]["x_full_df"].index
    )
    learning_xyz[condition]["y_full"][subround_id] = y_full

    learning_xyz[condition]["y_full_up"][subround_id] = y_full_up

    learning_xyz[condition]["y_train_df"][subround_id] = pd.DataFrame(
        y_train, index=learning_xyz[condition]["x_train_df"].index
    )
    learning_xyz[condition]["y_train"][subround_id] = y_train

    learning_xyz[condition]["y_train_up"][subround_id] = y_train_up

    learning_xyz[condition]["y_train_mixed_up"][subround_id] = y_train_mixed_up

    learning_xyz[condition]["y_test"][subround_id] = y_test

    return learning_xyz


def create_learninglist(
    learning_xyz: dict[str, dict[str, dict]],
    fract_full: dict[str, pd.DataFrame],
    fract_full_up: dict[str, pd.DataFrame],
    fract_marker: dict[str, pd.DataFrame],
    fract_marker_up: dict[str, pd.DataFrame],
    fract_test: dict[str, pd.DataFrame],
    condition: str,
    roundn: int,
):
    round_id = f"ROUND_{roundn}"
    classes = fract_marker[condition]["class"].unique().tolist()
    learning_xyz[condition]["classes"] = classes

    # TODO(performance): not all of those need to be stored or be converted
    #  to lists
    learning_xyz[condition]["W_full_df"] = fract_full[condition]["class"]
    learning_xyz[condition]["W_full"] = list(
        learning_xyz[condition]["W_full_df"]
    )
    learning_xyz[condition]["W_full_up_df"][round_id] = fract_full_up[
        condition
    ]["class"]
    learning_xyz[condition]["W_full_up"][round_id] = list(
        learning_xyz[condition]["W_full_up_df"][round_id]
    )
    learning_xyz[condition]["W_train_df"] = fract_marker[condition]["class"]
    learning_xyz[condition]["W_train"] = list(
        learning_xyz[condition]["W_train_df"]
    )
    learning_xyz[condition]["W_train_up_df"][round_id] = fract_marker_up[
        condition
    ]["class"]
    learning_xyz[condition]["W_train_up"][round_id] = list(
        learning_xyz[condition]["W_train_up_df"][round_id]
    )

    learning_xyz[condition]["x_full_df"] = fract_full[condition].drop(
        columns=["class"]
    )
    learning_xyz[condition]["x_full"] = learning_xyz[condition][
        "x_full_df"
    ].to_numpy(dtype=float)
    learning_xyz[condition]["x_full_up_df"][round_id] = fract_full_up[
        condition
    ].drop(columns=["class"])
    learning_xyz[condition]["x_full_up"][round_id] = learning_xyz[condition][
        "x_full_up_df"
    ][round_id].to_numpy(dtype=float)
    learning_xyz[condition]["x_train_df"] = fract_marker[condition].drop(
        columns=["class"]
    )
    learning_xyz[condition]["x_train"] = learning_xyz[condition][
        "x_train_df"
    ].to_numpy(dtype=float)
    learning_xyz[condition]["x_train_up_df"][round_id] = fract_marker_up[
        condition
    ].drop(columns=["class"])
    learning_xyz[condition]["x_train_up"][round_id] = learning_xyz[condition][
        "x_train_up_df"
    ][round_id].to_numpy(dtype=float)
    learning_xyz[condition]["x_test_df"] = fract_test[condition].drop(
        columns=["class"]
    )
    learning_xyz[condition]["x_test"] = learning_xyz[condition][
        "x_test_df"
    ].to_numpy(dtype=float)

    learning_xyz[condition]["Z_train_df"] = pd.get_dummies(
        fract_marker[condition]["class"]
    )[learning_xyz[condition]["classes"]]
    learning_xyz[condition]["Z_train"] = learning_xyz[condition][
        "Z_train_df"
    ].to_numpy(dtype=float)
    learning_xyz[condition]["V_full_up"][round_id] = learning_xyz[condition][
        "x_full_up"
    ][round_id]

    return learning_xyz, classes


def mix_profiles(
    mix_steps: list[float],
    NN_params: NeuralNetworkParametersModel,
    fract_marker_up,
    fract_unmixed_up,
    fract_mixed_up,
    condition: str,
):
    """Create mixed profiles.

    Create pairwise combinations of profiles and mix them according to the
    `mix_steps`. Return a random sample of the mixed profiles.
    """
    class_list = list(set(list(fract_marker_up[condition]["class"])))
    combinations = [
        (a, b)
        for idx, a in enumerate(class_list)
        for b in class_list[idx + 1 :]
    ]

    fract_mixed_up[condition] = copy.deepcopy(fract_unmixed_up[condition])

    cur = 1
    for comb in combinations:
        profiles_own = (
            fract_marker_up[condition]
            .copy()
            .loc[fract_marker_up[condition]["class"] == comb[0]]
            .drop(columns=["class"])
        )
        profiles_other = (
            fract_marker_up[condition]
            .copy()
            .loc[fract_marker_up[condition]["class"] == comb[1]]
            .drop(columns=["class"])
        )

        new_index = [
            f"{i}_{j}"
            for i, j in zip(profiles_own.index, profiles_other.index)
        ]
        for part in mix_steps:
            new_index_part = [
                f"{i + cur}_{value}" for i, value in enumerate(new_index)
            ]
            own_part = profiles_own.multiply(part)
            other_part = profiles_other.multiply(1 - part)

            own_part.index = new_index_part
            other_part.index = new_index_part

            profiles_mixed = own_part + other_part
            for classname in class_list:
                if classname == comb[0]:
                    profiles_mixed[classname] = part
                elif classname == comb[1]:
                    profiles_mixed[classname] = 1 - part
                else:
                    profiles_mixed[classname] = 0.0
            profiles_mixed = profiles_mixed.sample(frac=NN_params.mixed_batch)
            fract_mixed_up[condition] = pd.concat(
                [fract_mixed_up[condition], profiles_mixed]
            ).sample(frac=1)
            cur += len(profiles_mixed)
    return fract_mixed_up


def single_prediction(
    learning_xyz,
    clf: svm.SVC,
    svm_metrics,
    fract_marker,
    svm_marker,
    fract_test,
    svm_test,
    condition,
    roundn: int,
):
    """Perform single prediction.

    Train Support Vector Machine (SVM) classifier and predict the classes.


    :param learning_xyz: The learning data. This will be updated in place.
    """
    logger.info(f"Performing single prediction for {condition}...")
    round_id = f"ROUND_{roundn}"
    x_full = learning_xyz["x_full"]
    x_train = learning_xyz["x_train"]
    x_train_up = learning_xyz["x_train_up"][round_id]
    x_test = learning_xyz["x_test"]

    W_train = learning_xyz["W_train"]
    W_train_up = learning_xyz["W_train_up"][round_id]

    clf.fit(x_train_up, W_train_up)

    w_full = clf.predict(x_full).tolist()
    w_train = clf.predict(x_train).tolist()
    w_test = clf.predict(x_test).tolist()

    w_full_prob = list(map(max, list(clf.predict_proba(x_full))))
    w_train_prob = list(map(max, list(clf.predict_proba(x_train))))
    w_test_prob = list(map(max, list(clf.predict_proba(x_test))))

    confusion = pd.DataFrame(
        confusion_matrix(W_train, w_train, labels=list(clf.classes_)),
        index=clf.classes_,
        columns=clf.classes_,
    )
    accuracy = accuracy_score(W_train, w_train)
    precision = precision_score(W_train, w_train, average="macro")
    recall = recall_score(W_train, w_train, average="macro")
    f1 = f1_score(W_train, w_train, average="macro")

    svm_marker[condition] = copy.deepcopy(fract_marker[condition])
    svm_marker[condition]["svm_prediction"] = w_train
    svm_marker[condition]["svm_probability"] = w_train_prob

    svm_test[condition] = copy.deepcopy(fract_test[condition])
    svm_test[condition]["svm_prediction"] = w_test
    svm_test[condition]["svm_probability"] = w_test_prob

    svm_metrics[condition] = {
        "confusion": confusion,
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }

    learning_xyz["w_full"][round_id] = w_full
    learning_xyz["w_full_prob"][round_id] = w_full_prob
    learning_xyz["w_full_prob_df"][round_id] = copy.deepcopy(
        learning_xyz["x_full_df"]
    )
    learning_xyz["w_full_prob_df"][round_id]["SVM_winner"] = w_full
    learning_xyz["w_full_prob_df"][round_id]["SVM_prob"] = w_full_prob

    learning_xyz["w_train"][round_id] = w_train
    return svm_metrics, svm_marker, svm_test
