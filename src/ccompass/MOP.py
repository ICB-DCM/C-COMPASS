"""Multiple organelle prediction."""

import copy
import logging
from datetime import datetime
from typing import Any, Literal

import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
)

from ._utils import get_ccmps_data_directory
from .core import NeuralNetworkParametersModel, XYZ_Model

logger = logging.getLogger(__package__)


def upsample_condition(
    stds: pd.DataFrame,
    fract_full: pd.DataFrame,
    fract_marker: pd.DataFrame,
    method: Literal["none", "noised", "average", "noisedaverage"],
    noise_stds: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Perform upsampling for the given condition.

    :param method: The upsampling method to use.
    :param noise_stds: The noise level for upsampling in standard deviations.
        Only used for the "noised" and "noisedaverage" methods.
    :return: The upsampled marker and full profiles
    """
    fract_full_up = fract_full
    fract_marker_up = fract_marker.copy()

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
            if method == "noised":
                sample = data_class.sample(n=1)
                id_rnd = sample.index[0]
                name_up = f"up_{k}_{id_rnd}"
                k += 1

                profile_rnd_flat = sample.values.flatten()
                std_rnd = stds.loc[[id_rnd]]
                std_rnd = std_rnd[~std_rnd.index.duplicated(keep="first")]
                std_rnd_flat = std_rnd.values.flatten()
                std_rnd_flat = np.tile(
                    std_rnd_flat,
                    int(profile_rnd_flat.size / std_rnd_flat.size),
                )

                nv = np.random.normal(
                    profile_rnd_flat,
                    noise_stds * std_rnd_flat,
                    size=sample.shape,
                )
                nv = np.where(nv > 1, 1, np.where(nv < 0, 0, nv))
                profile_up = pd.DataFrame(nv, columns=sample.columns)

            elif method == "average":
                sample = data_class.sample(n=3, replace=True)
                name_up = f"up_{k}_{'_'.join(sample.index)}"
                k += 1
                profile_up = sample.median(axis=0).to_frame().transpose()

            elif method == "noisedaverage":
                sample = data_class.sample(n=3, replace=True)
                name_up = f"up_{k}_{'_'.join(sample.index)}"
                k += 1

                profile_av = sample.median(axis=0).to_frame().transpose()
                profile_av_flat = profile_av.values.flatten()
                nv = np.random.normal(
                    profile_av_flat,
                    noise_stds * class_std_flat,
                    size=profile_av.shape,
                )
                nv = np.where(nv > 1, 1, np.where(nv < 0, 0, nv))
                profile_up = pd.DataFrame(nv, columns=profile_av.columns)
            else:
                raise ValueError(f"Unknown upsampling method: {method}")

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


def MOP_exec(
    fract_full: dict[str, pd.DataFrame],
    fract_marker: dict[str, pd.DataFrame],
    fract_test: dict[str, pd.DataFrame],
    stds: dict[str, pd.DataFrame],
    nn_params: NeuralNetworkParametersModel,
) -> dict[str, XYZ_Model]:
    """Perform multi-organelle prediction.

    :param fract_full: dictionary of full profiles
    """
    conditions = list(fract_full.keys())
    learning_xyz = {condition: XYZ_Model() for condition in conditions}
    for condition in conditions:
        update_learninglist_const(
            learning_xyz[condition],
            fract_full[condition],
            fract_marker[condition],
            fract_test[condition],
        )

    # TODO: bring into shape for parallelization across conditions and rounds
    #  conditions are independent; rounds are independent
    #    only the rounds dict needs to be shared -> refactor to round_xyz model per condition
    #  for condition in conditions:
    for i_round in range(1, nn_params.rounds + 1):
        logger.info(f"Executing round {i_round}/{nn_params.rounds}...")
        round_id = f"ROUND_{i_round}"

        fract_full_up = {}
        fract_marker_up = {}

        if nn_params.upsampling:
            # upsample fractionation data for each condition x replicate
            for condition in conditions:
                logger.info(f"Upsampling condition {condition}")
                fract_marker_up[condition], fract_full_up[condition] = (
                    upsample_condition(
                        stds.get(condition),
                        fract_full[condition],
                        fract_marker[condition],
                        method=nn_params.upsampling_method,
                        noise_stds=nn_params.upsampling_noise,
                    )
                )
            logger.info("upsampling done!")
        else:
            fract_marker_up = copy.deepcopy(fract_marker)
            fract_full_up = copy.deepcopy(fract_full)

        for condition in conditions:
            update_learninglist_round(
                learning_xyz[condition],
                fract_full_up[condition],
                fract_marker_up[condition],
                round_id,
            )

        svm_marker = {}
        for condition in conditions:
            logger.info(f"Performing single prediction for {condition}...")
            _, svm_marker[condition], _ = single_prediction(
                learning_xyz[condition],
                fract_marker[condition],
                fract_test[condition],
                round_id,
            )

        if nn_params.svm_filter:
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
                        stds.get(condition),
                        fract_full[condition],
                        fract_marker_filtered[condition],
                        method=nn_params.upsampling_method,
                        noise_stds=nn_params.upsampling_noise,
                    )
                )
            logger.info("SVM filtering done.")

        fract_unmixed_up = {}
        for condition in conditions:
            # One-hot encode the classes
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

        if nn_params.mixed_part == "none":
            fract_mixed_up = copy.deepcopy(fract_unmixed_up)
        else:
            logger.info("Mixing profiles...")
            fract_mixed_up = {}
            mix_steps = [
                i / nn_params.mixed_part
                for i in range(1, nn_params.mixed_part)
            ]
            for condition in conditions:
                fract_mixed_up[condition] = mix_profiles(
                    fract_marker_up[condition],
                    fract_unmixed_up[condition],
                    mix_steps,
                    nn_params.mixed_batch,
                )

        for condition, xyz in learning_xyz.items():
            xyz.x_train_mixed_up_df[round_id] = fract_mixed_up[condition].drop(
                columns=xyz.classes
            )
            xyz.x_train_mixed_up[round_id] = xyz.x_train_mixed_up_df[
                round_id
            ].to_numpy(dtype=float)
            xyz.Z_train_mixed_up_df[round_id] = fract_mixed_up[condition][
                xyz.classes
            ]
            xyz.Z_train_mixed_up[round_id] = xyz.Z_train_mixed_up_df[
                round_id
            ].to_numpy(dtype=float)
        logger.info("mixing done!")

        for condition, xyz in learning_xyz.items():
            xyz.V_full_up[round_id] = xyz.x_full_up[round_id]

            if nn_params.AE == "none":
                # TODO(performance): Is there anything happening here?
                #  Just needless copying?
                y_full = xyz.x_full
                y_full_up = xyz.x_full_up[round_id]
                y_train = xyz.x_train
                y_train_up = xyz.x_train_up[round_id]
                y_train_mixed_up = xyz.x_train_mixed_up[round_id]
                y_test = xyz.x_test

                add_Y(
                    xyz,
                    y_full,
                    y_full_up,
                    y_train,
                    y_train_up,
                    y_train_mixed_up,
                    y_test,
                    i_round,
                    0,
                )
                for i_subround in range(1, nn_params.subrounds + 1):
                    add_Y(
                        xyz,
                        y_full,
                        y_full_up,
                        y_train,
                        y_train_up,
                        y_train_mixed_up,
                        y_test,
                        i_round,
                        i_subround,
                    )
            else:
                # TODO ADD AUTOENCODER HERE
                raise NotImplementedError("Autoencoder not implemented yet.")

        for condition in conditions:
            multi_predictions(
                learning_xyz[condition],
                nn_params,
                condition,
                i_round,
            )
    return learning_xyz


def multi_predictions(
    learning_xyz: XYZ_Model,
    nn_params: NeuralNetworkParametersModel,
    condition: str,
    roundn: int,
):
    """Perform multi organelle predictions."""
    logger.info(f"Training classifier for condition {condition}...")

    import keras_tuner as kt
    import tensorflow as tf

    from .classification_model import FNN_Classifier

    subround_id = f"ROUND_{roundn}_0"
    y_full = learning_xyz.y_full[subround_id]
    y_train = learning_xyz.y_train[subround_id]
    y_train_mixed_up = learning_xyz.y_train_mixed_up[subround_id]

    Z_train_mixed_up = learning_xyz.Z_train_mixed_up[f"ROUND_{roundn}"]
    set_shapes = [np.shape(y_train_mixed_up)[1], np.shape(Z_train_mixed_up)[1]]

    # Tune the hyperparameters
    classifier_directory = get_ccmps_data_directory()
    classifier_directory.mkdir(exist_ok=True, parents=True)

    now = datetime.now()
    time = now.strftime("%Y%m%d%H%M%S")
    tuner = kt.Hyperband(
        hypermodel=FNN_Classifier(set_shapes=set_shapes, nn_params=nn_params),
        hyperparameters=kt.HyperParameters(),
        objective="val_mean_squared_error",
        max_epochs=nn_params.NN_epochs,
        factor=3,
        directory=str(classifier_directory),
        project_name=f"{time}_Classifier_{condition}_{roundn}",
    )

    # noinspection PyUnresolvedReferences
    stop_early = tf.keras.callbacks.EarlyStopping(
        monitor="val_loss", patience=5
    )

    tuner.search(
        y_train_mixed_up,
        Z_train_mixed_up,
        epochs=nn_params.NN_epochs,
        validation_split=0.2,
        callbacks=[stop_early],
    )
    logger.info("Hyperparameter tuning done!")
    best_model = tuner.get_best_models(num_models=1)[0]
    best_hp = tuner.get_best_hyperparameters(num_trials=1)[0].values

    stringlist = []
    best_model.summary(print_fn=lambda x: stringlist.append(x))
    FNN_summary = "\n".join(stringlist)

    learning_xyz.FNN_summary[f"ROUND_{roundn}"] = FNN_summary

    z_full = best_model.predict(y_full)
    z_train = best_model.predict(y_train)

    add_Z(
        learning_xyz,
        z_full,
        z_train,
        roundn,
        0,
    )

    for subround in range(1, nn_params.subrounds + 1):
        logger.info(
            f"Training classifier for condition {condition} "
            f"{subround}/{nn_params.subrounds}..."
        )
        subround_id = f"ROUND_{roundn}_{subround}"
        y_full = learning_xyz.y_full[subround_id]
        y_train = learning_xyz.y_train[subround_id]
        y_train_mixed_up = learning_xyz.y_train_mixed_up[subround_id]

        fixed_model = FNN_Classifier(
            fixed_hp=best_hp, set_shapes=set_shapes, nn_params=nn_params
        ).build(None)
        fixed_model.fit(
            y_train_mixed_up,
            Z_train_mixed_up,
            epochs=nn_params.NN_epochs,
            validation_split=0.2,
            callbacks=[stop_early],
        )

        z_full = fixed_model.predict(y_full)
        z_train = fixed_model.predict(y_train)

        add_Z(
            learning_xyz,
            z_full,
            z_train,
            roundn,
            subround,
        )


def add_Z(
    learning_xyz: XYZ_Model,
    z_full,
    z_train,
    roundn: int,
    subroundn: int,
):
    subround_id = f"ROUND_{roundn}_{subroundn}"

    learning_xyz.z_full_df[subround_id] = pd.DataFrame(
        z_full,
        index=learning_xyz.y_full_df[subround_id].index,
        columns=learning_xyz.classes,
    )
    learning_xyz.z_full[subround_id] = z_full

    learning_xyz.z_train_df[subround_id] = pd.DataFrame(
        z_train,
        index=learning_xyz.y_train_df[subround_id].index,
        columns=learning_xyz.classes,
    )


def add_Y(
    learning_xyz: XYZ_Model,
    y_full,
    y_full_up,
    y_train,
    y_train_up,
    y_train_mixed_up,
    y_test,
    roundn: int,
    subroundn: int,
):
    subround_id = f"ROUND_{roundn}_{subroundn}"
    learning_xyz.y_full_df[subround_id] = pd.DataFrame(
        y_full, index=learning_xyz.x_full_df.index
    )
    learning_xyz.y_full[subround_id] = y_full
    learning_xyz.y_full_up[subround_id] = y_full_up
    learning_xyz.y_train_df[subround_id] = pd.DataFrame(
        y_train, index=learning_xyz.x_train_df.index
    )
    learning_xyz.y_train[subround_id] = y_train
    learning_xyz.y_train_up[subround_id] = y_train_up
    learning_xyz.y_train_mixed_up[subround_id] = y_train_mixed_up
    learning_xyz.y_test[subround_id] = y_test


def update_learninglist_const(
    learning_xyz: XYZ_Model,
    fract_full: pd.DataFrame,
    fract_marker: pd.DataFrame,
    fract_test: pd.DataFrame,
) -> None:
    """Populate `learning_xyz` with the learning data that is constant across
    rounds."""
    # TODO(performance): not all of those need to be stored or be converted
    #  to lists
    learning_xyz.classes = fract_marker["class"].unique().tolist()
    learning_xyz.W_full_df = fract_full["class"]
    learning_xyz.W_train_df = fract_marker["class"]
    learning_xyz.W_train = list(learning_xyz.W_train_df)
    learning_xyz.x_full_df = fract_full.drop(columns=["class"])
    learning_xyz.x_full = learning_xyz.x_full_df.to_numpy(dtype=float)
    learning_xyz.x_train_df = fract_marker.drop(columns=["class"])
    learning_xyz.x_train = learning_xyz.x_train_df.to_numpy(dtype=float)
    learning_xyz.x_test_df = fract_test.drop(columns=["class"])
    learning_xyz.x_test = learning_xyz.x_test_df.to_numpy(dtype=float)
    learning_xyz.Z_train_df = pd.get_dummies(fract_marker["class"])[
        learning_xyz.classes
    ]


def update_learninglist_round(
    learning_xyz: XYZ_Model,
    fract_full_up: pd.DataFrame,
    fract_marker_up: pd.DataFrame,
    round_id: str,
) -> None:
    """Populate `learning_xyz` with the learning data that is specific to the
    given round."""
    learning_xyz.W_train_up_df[round_id] = fract_marker_up["class"]
    learning_xyz.W_train_up[round_id] = list(
        learning_xyz.W_train_up_df[round_id]
    )

    learning_xyz.x_full_up_df[round_id] = fract_full_up.drop(columns=["class"])
    learning_xyz.x_full_up[round_id] = learning_xyz.x_full_up_df[
        round_id
    ].to_numpy(dtype=float)
    learning_xyz.x_train_up_df[round_id] = fract_marker_up.drop(
        columns=["class"]
    )
    learning_xyz.x_train_up[round_id] = learning_xyz.x_train_up_df[
        round_id
    ].to_numpy(dtype=float)
    learning_xyz.V_full_up[round_id] = learning_xyz.x_full_up[round_id]


def mix_profiles(
    fract_marker_up,
    fract_unmixed_up,
    mix_steps: list[float],
    mixed_batch: float,
) -> pd.DataFrame:
    """Create mixed profiles.

    Create pairwise combinations of profiles and mix them according to the
    `mix_steps`. Return a random sample of the mixed profiles.
    """
    class_list = list(set(list(fract_marker_up["class"])))
    combinations = [
        (a, b)
        for idx, a in enumerate(class_list)
        for b in class_list[idx + 1 :]
    ]

    fract_mixed_up = copy.deepcopy(fract_unmixed_up)

    cur = 1
    for comb in combinations:
        # TODO (performance): We can avoid some copying here
        profiles_own = (
            fract_marker_up.copy()
            .loc[fract_marker_up["class"] == comb[0]]
            .drop(columns=["class"])
        )
        profiles_other = (
            fract_marker_up.copy()
            .loc[fract_marker_up["class"] == comb[1]]
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
            profiles_mixed = profiles_mixed.sample(frac=mixed_batch)
            fract_mixed_up = pd.concat(
                [fract_mixed_up, profiles_mixed]
            ).sample(frac=1)
            cur += len(profiles_mixed)
    return fract_mixed_up


def single_prediction(
    learning_xyz: XYZ_Model,
    fract_marker: pd.DataFrame,
    fract_test: pd.DataFrame,
    round_id: str,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """Perform single-class (single-compartment) prediction.

    Train Support Vector Machine (SVM) classifier and predict the classes.

    :param learning_xyz: The learning data. This will be updated in place.
    """
    x_full = learning_xyz.x_full
    x_train = learning_xyz.x_train
    x_train_up = learning_xyz.x_train_up[round_id]
    x_test = learning_xyz.x_test
    W_train = learning_xyz.W_train
    W_train_up = learning_xyz.W_train_up[round_id]

    # train classifier on the upsampled data
    clf = svm.SVC(kernel="rbf", probability=True)
    clf.fit(x_train_up, W_train_up)

    # predict the classes
    # TODO(performance): No need to predict x_full,
    #  since that is x_train + x_test
    w_full = clf.predict(x_full).tolist()
    w_train = clf.predict(x_train).tolist()
    w_test = clf.predict(x_test).tolist()

    w_full_prob = clf.predict_proba(x_full).max(axis=1)
    w_train_prob = clf.predict_proba(x_train).max(axis=1)
    w_test_prob = clf.predict_proba(x_test).max(axis=1)

    confusion = pd.DataFrame(
        confusion_matrix(W_train, w_train, labels=list(clf.classes_)),
        index=clf.classes_,
        columns=clf.classes_,
    )
    accuracy = accuracy_score(W_train, w_train)
    precision = precision_score(W_train, w_train, average="macro")
    recall = recall_score(W_train, w_train, average="macro")
    f1 = f1_score(W_train, w_train, average="macro")

    svm_marker = copy.deepcopy(fract_marker)
    svm_marker["svm_prediction"] = w_train
    svm_marker["svm_probability"] = w_train_prob

    svm_test = copy.deepcopy(fract_test)
    svm_test["svm_prediction"] = w_test
    svm_test["svm_probability"] = w_test_prob

    svm_metrics = {
        "confusion": confusion,
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }

    learning_xyz.w_full[round_id] = w_full
    learning_xyz.w_full_prob[round_id] = w_full_prob
    learning_xyz.w_full_prob_df[round_id] = copy.deepcopy(
        learning_xyz.x_full_df
    )
    learning_xyz.w_full_prob_df[round_id]["SVM_winner"] = w_full
    learning_xyz.w_full_prob_df[round_id]["SVM_prob"] = w_full_prob
    learning_xyz.w_train[round_id] = w_train

    return svm_metrics, svm_marker, svm_test
