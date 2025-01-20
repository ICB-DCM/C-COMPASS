"""Multiple organelle prediction."""

import copy
import random
from datetime import datetime

import keras.backend as K
import keras_tuner as kt
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow.keras
from keras import ops
from scipy import stats
from sklearn import svm
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
)
from tensorflow import keras

from ._utils import get_ccmps_data_directory
from .core import NeuralNetworkParametersModel

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
            # units_init = np.shape(y_train_mixed_up)[1]
            # model.add(tf.keras.Input(units_init,))

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


def upsampling(
    NN_params: NeuralNetworkParametersModel,
    stds,
    fract_full,
    fract_full_up,
    fract_marker,
    fract_marker_up,
    condition,
):
    """Perform upsampling for all conditions."""
    fract_full_up[condition] = fract_full[condition]
    fract_marker_up[condition] = fract_marker[condition]

    if NN_params.upsampling_method == "none":
        pass
    else:
        print(f"Upsampling condition {condition}:")
        class_sizes = {}
        for classname in list(set(fract_marker[condition]["class"])):
            class_sizes[classname] = list(
                fract_marker[condition]["class"]
            ).count(classname)
        class_maxsize = max(class_sizes.values())
        k = 1
        for classname in list(set(fract_marker[condition]["class"])):
            print(classname)
            data_class_temp = fract_marker[condition].loc[
                fract_marker[condition]["class"] == classname
            ]
            data_class = data_class_temp.drop(columns=["class"])
            class_difference = abs(class_maxsize - class_sizes[classname])

            if class_sizes[classname] > class_maxsize:
                ID_rnd = random.sample(
                    list(data_class.index), class_difference - 1
                )
                fract_marker_up[condition].drop(ID_rnd, inplace=True)
            if class_sizes[classname] < class_maxsize:
                class_up = pd.DataFrame(columns=data_class.columns)

                class_std = data_class.std(axis=0).to_frame().transpose()
                class_std_flat = class_std.values.flatten()

                for i in range(class_difference):
                    if NN_params.upsampling_method == "noised":
                        ID_rnd = random.choice(list(data_class.index))
                        name_up = f"up_{k}_{ID_rnd}"
                        k += 1

                        profile_rnd = data_class.loc[[ID_rnd]]
                        profile_rnd = profile_rnd[
                            ~profile_rnd.index.duplicated(keep="first")
                        ]
                        profile_rnd_flat = profile_rnd.values.flatten()

                        std_rnd = stds[condition].loc[[ID_rnd]]
                        std_rnd = std_rnd[
                            ~std_rnd.index.duplicated(keep="first")
                        ]
                        std_rnd_flat = std_rnd.values.flatten()
                        std_rnd_flat = np.tile(
                            std_rnd_flat,
                            int(profile_rnd_flat.size / std_rnd_flat.size),
                        )

                        nv = np.random.normal(
                            profile_rnd_flat,
                            NN_params.upsampling_noise * std_rnd_flat,
                            size=profile_rnd.shape,
                        )
                        nv = np.where(nv > 1, 1, np.where(nv < 0, 0, nv))

                        profile_up = pd.DataFrame(
                            nv, columns=profile_rnd.columns
                        )

                    elif NN_params.upsampling_method == "average":
                        ID_rnd_1 = random.choice(list(data_class.index))
                        ID_rnd_2 = random.choice(list(data_class.index))
                        ID_rnd_3 = random.choice(list(data_class.index))
                        name_up = f"up_{k}_{ID_rnd_1}_{ID_rnd_2}_{ID_rnd_3}"
                        k += 1

                        profile_rnd_1 = data_class.loc[[ID_rnd_1]]
                        profile_rnd_1 = profile_rnd_1[
                            ~profile_rnd_1.index.duplicated(keep="first")
                        ]
                        profile_rnd_2 = data_class.loc[[ID_rnd_2]]
                        profile_rnd_2 = profile_rnd_2[
                            ~profile_rnd_2.index.duplicated(keep="first")
                        ]
                        profile_rnd_3 = data_class.loc[[ID_rnd_3]]
                        profile_rnd_3 = profile_rnd_3[
                            ~profile_rnd_3.index.duplicated(keep="first")
                        ]

                        profile_up = (
                            pd.concat(
                                [profile_rnd_1, profile_rnd_2, profile_rnd_3]
                            )
                            .median(axis=0)
                            .to_frame()
                            .transpose()
                        )

                    elif NN_params.upsampling_method == "noisedaverage":
                        ID_rnd_1 = random.choice(list(data_class.index))
                        ID_rnd_2 = random.choice(list(data_class.index))
                        ID_rnd_3 = random.choice(list(data_class.index))
                        name_up = f"up_{k}_{ID_rnd_1}_{ID_rnd_2}_{ID_rnd_3}"
                        k += 1

                        profile_rnd_1 = data_class.loc[[ID_rnd_1]]
                        profile_rnd_1 = profile_rnd_1[
                            ~profile_rnd_1.index.duplicated(keep="first")
                        ]
                        profile_rnd_2 = data_class.loc[[ID_rnd_2]]
                        profile_rnd_2 = profile_rnd_2[
                            ~profile_rnd_2.index.duplicated(keep="first")
                        ]
                        profile_rnd_3 = data_class.loc[[ID_rnd_3]]
                        profile_rnd_3 = profile_rnd_3[
                            ~profile_rnd_3.index.duplicated(keep="first")
                        ]

                        profile_av = (
                            pd.concat(
                                [profile_rnd_1, profile_rnd_2, profile_rnd_3]
                            )
                            .median(axis=0)
                            .to_frame()
                            .transpose()
                        )
                        profile_av_flat = profile_av.values.flatten()

                        nv = np.random.normal(
                            profile_av_flat,
                            NN_params.upsampling_noise * class_std_flat,
                            size=profile_av.shape,
                        )
                        nv = np.where(nv > 1, 1, np.where(nv < 0, 0, nv))

                        profile_up = pd.DataFrame(
                            nv, columns=profile_av.columns
                        )
                    else:
                        raise ValueError(
                            f"Unknown upsampling method: {NN_params.upsampling_method}"
                        )

                    profile_up.index = [name_up]
                    profile_up["class"] = [classname]

                    class_up = pd.concat(
                        [class_up, profile_up], axis=0, ignore_index=False
                    )

                fract_marker_up[condition] = pd.concat(
                    [fract_marker_up[condition], class_up],
                    axis=0,
                    ignore_index=False,
                )
                # print(len(class_up))
                # print(len(fract_marker_up[condition]['class']))
                # print(fract_marker_up[condition])
                fract_full_up[condition] = pd.concat(
                    [fract_full_up[condition], class_up],
                    axis=0,
                    ignore_index=False,
                )
        print("")
    fract_marker_up[condition] = fract_marker_up[condition].sample(frac=1)
    fract_full_up[condition] = fract_full_up[condition].sample(frac=1)
    return fract_marker_up, fract_full_up


def create_fullprofiles(
    fract_marker: dict[str, pd.DataFrame], fract_test: dict[str, pd.DataFrame]
) -> dict[str, pd.DataFrame]:
    fract_full = {}
    for condition in fract_test:
        fract_full[condition] = pd.concat(
            [fract_test[condition], fract_marker[condition]]
        )
    return fract_full


def sum1_normalization(x):
    """Normalize the input to sum to 1."""
    return x / (ops.sum(x, axis=1, keepdims=True) + K.epsilon())


def combine_rounds(data):
    dfs = []
    for roundn in data:
        dfs.append(data[roundn])

    concatenated_df = pd.concat(dfs, axis=0)

    mean_df = concatenated_df.groupby(level=0).mean()
    std_df = concatenated_df.groupby(level=0).std()

    count_non_na = concatenated_df.groupby(level=0).apply(
        lambda x: x.notna().sum()
    )

    count_df = pd.DataFrame(index=mean_df.index)
    count_df["count"] = count_non_na.max(axis=1)

    return mean_df, std_df, count_df


def pairwise_t_test(
    conditions_values, conditions_std, conditions_count, rounds, alpha=0.05
):
    significant_changes = []
    num_conditions = len(conditions_values)

    for i in range(num_conditions):
        for j in range(i + 1, num_conditions):
            values_df1, values_df2 = conditions_values[i], conditions_values[j]
            std_df1, std_df2 = conditions_std[i], conditions_std[j]
            count_df1, count_df2 = conditions_count[i], conditions_count[j]

            common_indices = values_df1.index.intersection(values_df2.index)

            for idx in common_indices:
                for col in values_df1.columns:
                    mean1, std1 = values_df1.at[idx, col], std_df1.at[idx, col]
                    mean2, std2 = values_df2.at[idx, col], std_df2.at[idx, col]
                    n1, n2 = (
                        count_df1.at[idx, "count"] / rounds,
                        count_df2.at[idx, "count"] / rounds,
                    )

                    t_stat, p_value = stats.ttest_ind_from_stats(
                        mean1, std1, n1, mean2, std2, n2
                    )

                    if p_value < alpha:
                        magnitude_of_change = mean1 - mean2
                        significant_changes.append(
                            (
                                idx,
                                col,
                                f"Condition {i + 1} vs {j + 1}",
                                t_stat,
                                p_value,
                                magnitude_of_change,
                            )
                        )

    return pd.DataFrame(
        significant_changes,
        columns=[
            "Index",
            "Class",
            "Comparison",
            "T-Stat",
            "P-Value",
            "MagnitudeOfChange",
        ],
    )


# ----------------------------------------------------------------------------


def MOP_exec(
    fract_conditions,
    fract_full,
    fract_marker_old,
    fract_test,
    fract_std,
    fract_info,
    key,
    NN_params: NeuralNetworkParametersModel,
):
    """Perform multi-organelle prediction."""
    conditions_std = [x for x in fract_conditions if x != "[KEEP]"]
    conditions = [x for x in fract_full]

    stds = {}
    if not key == "[IDENTIFIER]":
        for condition in conditions_std:
            stds[condition] = pd.merge(
                fract_std["class"][condition],
                fract_info[key],
                left_index=True,
                right_index=True,
                how="left",
            ).set_index(key)

    # -------------------------
    ## UPSAMPLING START
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
        # learning_xyz[condition]['Z_train_up_df'] = {}
        learning_xyz[condition]["Z_train_up"] = {}
        # learning_xyz[condition]['V_full_up_df'] = {}
        learning_xyz[condition]["V_full_up"] = {}
        learning_xyz[condition]["x_train_mixed_up_df"] = {}
        learning_xyz[condition]["x_train_mixed_up"] = {}
        learning_xyz[condition]["Z_train_mixed_up_df"] = {}
        learning_xyz[condition]["Z_train_mixed_up"] = {}
        learning_xyz[condition]["AE_summary"] = {}
        learning_xyz[condition]["AE_history"] = {}
        # learning_xyz[condition]['v_full_df'] = {}
        # learning_xyz[condition]['v_full'] = {}
        # learning_xyz[condition]['v_full_up_df'] = {}
        # learning_xyz[condition]['v_full_up'] = {}
        learning_xyz[condition]["y_full_df"] = {}
        learning_xyz[condition]["y_full"] = {}
        # learning_xyz[condition]['y_full_up_df'] = {}
        learning_xyz[condition]["y_full_up"] = {}
        learning_xyz[condition]["y_train_df"] = {}
        learning_xyz[condition]["y_train"] = {}
        # learning_xyz[condition]['y_train_up_df'] = {}
        learning_xyz[condition]["y_train_up"] = {}
        # learning_xyz[condition]['y_train_mixed_up_df'] = {}
        learning_xyz[condition]["y_train_mixed_up"] = {}
        # learning_xyz[condition]['y_test_df'] = {}
        learning_xyz[condition]["y_test"] = {}
        learning_xyz[condition]["FNN_summary"] = {}
        learning_xyz[condition]["FNN_history"] = {}
        learning_xyz[condition]["z_full_df"] = {}
        learning_xyz[condition]["z_full"] = {}
        learning_xyz[condition]["z_train_df"] = {}
        learning_xyz[condition]["z_train"] = {}
        # learning_xyz[condition]['z_train_up_df'] = {}
        # learning_xyz[condition]['z_train_up'] = {}
        # learning_xyz[condition]['z_train_mixed_up_df'] = {}
        # learning_xyz[condition]['z_train_mixed_up'] = {}
        # learning_xyz[condition]['z_test_df'] = {}
        # learning_xyz[condition]['z_test'] = {}

    for R in range(1, NN_params.rounds + 1):
        # print('ROUND 0')
        print("upsampling...")
        # fract_full_up = copy.deepcopy(fract_full)
        # fract_marker_up = copy.deepcopy(fract_marker_old)

        fract_full_up = {}
        fract_marker_up = {}

        fract_marker = copy.deepcopy(fract_marker_old)

        if NN_params.upsampling:
            for condition in conditions:
                fract_marker_up, fract_full_up = upsampling(
                    NN_params,
                    stds,
                    fract_full,
                    fract_full_up,
                    fract_marker,
                    fract_marker_up,
                    condition,
                )
        else:
            fract_marker_up = copy.deepcopy(fract_marker)
            fract_full_up = copy.deepcopy(fract_full)

        print("upsampling done!")
        print("")

        print("creating data...")
        for condition in conditions:
            learning_xyz, classes = create_learninglist(
                learning_xyz,
                fract_full,
                fract_full_up,
                fract_marker,
                fract_marker_up,
                fract_test,
                condition,
                R,
                0,
            )

        print("data complete!")
        print("")

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
                R,
            )

        if NN_params.svm_filter:
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
                # fract_full_up = copy.deepcopy(fract_full)
                # fract_marker_up = copy.deepcopy(fract_marker_old)
                # fract_marker = copy.deepcopy(fract_marker_old)
                fract_marker_up, fract_full_up = upsampling(
                    NN_params,
                    stds,
                    fract_full,
                    fract_full_up,
                    fract_marker_filtered,
                    fract_marker_up,
                    condition,
                )

        fract_unmixed_up = {}
        for condition in conditions:
            # print()
            unmixed_dummies = pd.get_dummies(
                fract_marker_up[condition]["class"]
            )  # [learning_xyz[condition]['classes']]
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
            fract_mixed_up = {}
            mix_steps = [
                i / (NN_params.mixed_part)
                for i in range(1, (NN_params.mixed_part))
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
        round_id = f"ROUND_{R}"
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
        print("mixing done!")
        print("")

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
                    R,
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
                        R,
                        SR,
                    )
            else:
                # TODO ADD AUTOENCODER HERE
                raise NotImplementedError("Autoencoder not implemented yet.")

        FNN_classifier = _create_classifier_hypermodel(NN_params)
        for condition in conditions:
            FNN_ens, learning_xyz = multi_predictions(
                FNN_classifier, learning_xyz, NN_params, condition, R
            )

    # print('post-processing...')
    # for condition in conditions:
    #     print(condition)
    #     learning_xyz[condition]['z_full_mean_df'], learning_xyz[condition]['z_full_std_df'], learning_xyz[condition]['n_full_df'] = combine_rounds(learning_xyz[condition]['z_full_df'])
    #     learning_xyz[condition]['z_train_mean_df'], learning_xyz[condition]['z_train_std_df'], learning_xyz[condition]['n_train_df'] = combine_rounds(learning_xyz[condition]['z_train_df'])
    #     learning_xyz[condition]['z_train_mixed_up_mean_df'], learning_xyz[condition]['z_train_mixed_up_std_df'], learning_xyz[condition]['n_train_mixed_up_df'] = combine_rounds(learning_xyz[condition]['z_train_mixed_up_df'])
    #     learning_xyz[condition]['z_test_mean_df'], learning_xyz[condition]['z_test_std_df'], learning_xyz[condition]['n_test_df'] = combine_rounds(learning_xyz[condition]['z_test_df'])

    #     learning_xyz[condition]['z_full_mean_filtered_df'] = learning_xyz[condition]['z_full_mean_df']
    #     learning_xyz[condition]['z_train_mean_filtered_df'] = learning_xyz[condition]['z_train_mean_df']
    #     learning_xyz[condition]['z_test_mean_filtered_df'] = learning_xyz[condition]['z_test_mean_df']
    #     learning_xyz[condition]['Z_train_df'] = learning_xyz[condition]['Z_train_df'][~learning_xyz[condition]['Z_train_df'].index.duplicated(keep='first')]
    #     for class_act in learning_xyz[condition]['classes']:
    #         nonmarker_pred = learning_xyz[condition]['z_train_mean_df'].loc[learning_xyz[condition]['Z_train_df'][class_act] != 1]
    #         thresh = np.percentile(nonmarker_pred[class_act].tolist(), NN_params['reliability'])
    #         learning_xyz[condition]['z_full_mean_filtered_df'].loc[learning_xyz[condition]['z_full_mean_filtered_df'][class_act] < thresh, class_act] = 0.
    #         learning_xyz[condition]['z_train_mean_filtered_df'].loc[learning_xyz[condition]['z_train_mean_filtered_df'][class_act] < thresh, class_act] = 0.
    #         learning_xyz[condition]['z_test_mean_filtered_df'].loc[learning_xyz[condition]['z_test_mean_filtered_df'][class_act] < thresh, class_act] = 0.

    # print('')
    # print('making statistics...')
    # prediction_train = []
    # prediction_test = []
    # std_train = []
    # std_test = []
    # count_train = []
    # count_test = []
    # for condition in conditions:
    #     print(condition)
    #     prediction_train.append(learning_xyz[condition]['z_train_mean_filtered_df'])
    #     prediction_test.append(learning_xyz[condition]['z_test_mean_filtered_df'])
    #     std_train.append(learning_xyz[condition]['z_train_std_df'])
    #     std_test.append(learning_xyz[condition]['z_test_std_df'])
    #     count_train.append(learning_xyz[condition]['n_train_df'])
    #     count_test.append(learning_xyz[condition]['n_test_df'])

    # print('testing significance...')
    # significant_train = pairwise_t_test(prediction_train, std_train, count_train, NN_params['rounds'])
    # significant_test = pairwise_t_test(prediction_test, std_test, count_test, NN_params['rounds'])
    # print('DONE!')

    # now = datetime.now()
    # print(now.strftime("%Y%m%d%H%M%S"))
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
    y_full = learning_xyz[condition]["y_full"][f"ROUND_{roundn}_0"]
    y_train = learning_xyz[condition]["y_train"][f"ROUND_{roundn}_0"]
    y_train_up = learning_xyz[condition]["y_train_up"][f"ROUND_{roundn}_0"]
    y_train_mixed_up = learning_xyz[condition]["y_train_mixed_up"][
        f"ROUND_{roundn}_0"
    ]
    y_test = learning_xyz[condition]["y_test"][f"ROUND_{roundn}_0"]

    Z_train_mixed_up = learning_xyz[condition]["Z_train_mixed_up"][
        f"ROUND_{roundn}"
    ]
    set_shapes = [np.shape(y_train_mixed_up)[1], np.shape(Z_train_mixed_up)[1]]

    # Tune the hyperparameters
    classifier_directory = get_ccmps_data_directory()
    classifier_directory.mkdir(exist_ok=True, parents=True)

    now = datetime.now()
    time = now.strftime("%Y%m%d%H%M%S")
    FNN_tuner = kt.Hyperband(
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

    FNN_tuner.search(
        y_train_mixed_up,
        Z_train_mixed_up,
        epochs=NN_params.NN_epochs,
        validation_split=0.2,
        callbacks=[stop_early],
    )

    FNN_best = FNN_tuner.get_best_models(num_models=1)[0]
    best_hp = FNN_tuner.get_best_hyperparameters(num_trials=1)[0].values

    FNN_ens = [FNN_best]
    # FNN_ens[0] = copy.deepcopy(FNN_best)
    # FNN_best.build(y_train_mixed_up.shape)

    stringlist = []
    FNN_best.summary(print_fn=lambda x: stringlist.append(x))
    FNN_summary = "\n".join(stringlist)

    learning_xyz[condition]["FNN_summary"]["ROUND_roundn"] = FNN_summary

    z_full = FNN_best.predict(y_full)
    z_train = FNN_best.predict(y_train)
    z_train_up = FNN_best.predict(y_train_up)
    z_train_mixed_up = FNN_best.predict(y_train_mixed_up)
    z_test = FNN_best.predict(y_test)

    learning_xyz = add_Z(
        learning_xyz,
        z_full,
        z_train,
        z_train_up,
        z_train_mixed_up,
        z_test,
        condition,
        roundn,
        0,
    )

    for subround in range(1, NN_params.subrounds + 1):
        # print(learning_xyz[condition]['y_full'])
        subround_id = f"ROUND_{roundn}_{subround}"
        y_full = learning_xyz[condition]["y_full"][subround_id]
        y_train = learning_xyz[condition]["y_train"][subround_id]
        y_train_up = learning_xyz[condition]["y_train_up"][subround_id]
        y_train_mixed_up = learning_xyz[condition]["y_train_mixed_up"][
            subround_id
        ]
        y_test = learning_xyz[condition]["y_test"][subround_id]

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
        FNN_ens.append(fixed_model)

        z_full = fixed_model.predict(y_full)
        z_train = fixed_model.predict(y_train)
        z_train_up = fixed_model.predict(y_train_up)
        z_train_mixed_up = fixed_model.predict(y_train_mixed_up)
        z_test = fixed_model.predict(y_test)

        learning_xyz = add_Z(
            learning_xyz,
            z_full,
            z_train,
            z_train_up,
            Z_train_mixed_up,
            z_test,
            condition,
            roundn,
            subround,
        )

    return FNN_ens, learning_xyz


def add_Z(
    learning_xyz,
    z_full,
    z_train,
    z_train_up,
    z_train_mixed_up,
    z_test,
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

    # learning_xyz[condition]['z_train_up_df']['ROUND_' + str(roundn) + '_' + str(subroundn)] = pd.DataFrame(z_train_up,
    #                                                                            index = learning_xyz[condition]['y_train_up_df']['ROUND_' + str(roundn) + '_' + str(subroundn)].index,
    #                                                                            columns = learning_xyz[condition]['classes'])
    # learning_xyz[condition]['z_train_up']['ROUND_' + str(roundn) + '_' + str(subroundn)] = z_train_up

    # learning_xyz[condition]['z_train_mixed_up_df']['ROUND_' + str(roundn) + '_' + str(subroundn)] = pd.DataFrame(z_train_mixed_up,
    #                                                                                  index = learning_xyz[condition]['y_train_mixed_up_df']['ROUND_' + str(roundn) + '_' + str(subroundn)].index,
    #                                                                                  columns = learning_xyz[condition]['classes'])
    # learning_xyz[condition]['z_train_mixed_up']['ROUND_' + str(roundn) + '_' + str(subroundn)] = z_train_mixed_up

    # learning_xyz[condition]['z_test_df']['ROUND_' + str(roundn) + '_' + str(subroundn)] = pd.DataFrame(z_test,
    #                                                                        index = learning_xyz[condition]['y_test_df']['ROUND_' + str(roundn) + '_' + str(subroundn)].index,
    #                                                                        columns = learning_xyz[condition]['classes'])
    # learning_xyz[condition]['z_test']['ROUND_' + str(roundn) + '_' + str(subroundn)] = z_test

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

    # learning_xyz[condition]['y_full_up_df']['ROUND_' + str(roundn) + '_' + str(subroundn)] = pd.DataFrame(y_full_up,
    #                                                                           index = learning_xyz[condition]['x_full_up_df']['ROUND_' + str(roundn)].index)
    learning_xyz[condition]["y_full_up"][subround_id] = y_full_up

    learning_xyz[condition]["y_train_df"][subround_id] = pd.DataFrame(
        y_train, index=learning_xyz[condition]["x_train_df"].index
    )
    learning_xyz[condition]["y_train"][subround_id] = y_train

    # learning_xyz[condition]['y_train_up_df']['ROUND_' + str(roundn) + '_' + str(subroundn)] = pd.DataFrame(y_train_up,
    #                                                                            index = learning_xyz[condition]['x_train_up_df']['ROUND_' + str(roundn)].index)
    learning_xyz[condition]["y_train_up"][subround_id] = y_train_up

    # learning_xyz[condition]['y_train_mixed_up_df']['ROUND_' + str(roundn) + '_' + str(subroundn)] = pd.DataFrame(y_train_mixed_up,
    #                                                                                  index = learning_xyz[condition]['x_train_mixed_up_df']['ROUND_' + str(roundn)].index)
    learning_xyz[condition]["y_train_mixed_up"][subround_id] = y_train_mixed_up

    # learning_xyz[condition]['y_test_df']['ROUND_' + str(roundn) + '_' + str(subroundn)] = pd.DataFrame(y_test,
    #                                                                        index = learning_xyz[condition]['x_test_df'].index)
    learning_xyz[condition]["y_test"][subround_id] = y_test

    return learning_xyz


def create_learninglist(
    learning_xyz,
    fract_full,
    fract_full_up,
    fract_marker,
    fract_marker_up,
    fract_test,
    condition,
    roundn: int,
    subroundn: int,
):
    round_id = f"ROUND_{roundn}"
    classes = list(set(fract_marker[condition]["class"]))
    learning_xyz[condition]["classes"] = classes

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
    # learning_xyz[condition]['Z_train_up_df']['ROUND_' + str(roundn) + '_' + str(subroundn)] = pd.get_dummies(fract_marker_up[condition]['class'])[learning_xyz[condition]['classes']]
    # learning_xyz[condition]['Z_train_up']['ROUND_' + str(roundn) + '_' + str(subroundn)] = learning_xyz[condition]['Z_train_up_df']['ROUND_' + str(roundn) + '_' + str(subroundn)].to_numpy(dtype = float)

    # learning_xyz[condition]['V_full_up_df']['ROUND_' + str(roundn)] = learning_xyz[condition]['x_full_up_df']['ROUND_' + str(roundn)]
    learning_xyz[condition]["V_full_up"][round_id] = learning_xyz[condition][
        "x_full_up"
    ][round_id]

    return learning_xyz, classes


def mix_profiles(
    mix_steps,
    NN_params: NeuralNetworkParametersModel,
    fract_marker_up,
    fract_unmixed_up,
    fract_mixed_up,
    condition,
):
    class_list = list(set(list(fract_marker_up[condition]["class"])))
    combinations = [
        (a, b)
        for idx, a in enumerate(class_list)
        for b in class_list[idx + 1 :]
    ]

    fract_mixed_up[condition] = copy.deepcopy(fract_unmixed_up[condition])
    # print(condition)
    # print(fract_unmixed_up[condition])
    # print(fract_mixed_up[condition])

    cur = 1
    # print(fract_marker_up[condition])
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

    :param learning_xyz: The learning data. This will be updated in place.
    """
    print(condition)
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
    # learning_xyz['w_train_prob']['ROUND_' + str(roundn)] = w_train_prob

    # learning_xyz['w_test']['ROUND_' + str(roundn)] = w_test
    # learning_xyz['w_test_prob']['ROUND_' + str(roundn)] = w_test_prob

    return svm_metrics, svm_marker, svm_test
