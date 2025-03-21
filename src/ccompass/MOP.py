"""Multiple organelle prediction."""

import copy
import gc
import logging
import multiprocessing as mp
import queue
import threading
from collections.abc import Callable
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

from ._utils import (
    PrefixFilter,
    get_ccmps_data_directory,
    get_mp_ctx,
    stdout_to_logger,
)
from .core import (
    ConditionPredictionModel,
    NeuralNetworkParametersModel,
    TrainingRoundModel,
)

logger = logging.getLogger(__package__)


def upsample_condition(
    fract_marker: pd.DataFrame,
    method: Literal["none", "average", "noisedaverage"],
    noise_stds: float,
) -> pd.DataFrame:
    """Perform upsampling for the given condition.

    :param method: The upsampling method to use.
    :param noise_stds: The noise level for upsampling in standard deviations.
        Only used for the "noised" and "noisedaverage" methods.
    :return: The upsampled marker and full profiles
    """
    if fract_marker.empty:
        raise ValueError("Empty marker profile")

    fract_marker_up = fract_marker.copy()

    class_sizes = fract_marker["class"].value_counts()
    class_maxsize = class_sizes.max()
    logger.debug(
        f"Upsampling to N_max = {class_maxsize} "
        f"(class sizes: {dict(class_sizes)})"
    )
    k = 1
    for classname, data_class in fract_marker.groupby("class"):
        if not (class_difference := class_maxsize - class_sizes[classname]):
            continue

        data_class = data_class.drop(columns=["class"])
        class_up = pd.DataFrame(columns=data_class.columns)

        if method == "noisedaverage":
            class_std = data_class.std(axis=0).to_frame().transpose()
            class_std_flat = class_std.values.flatten()
        else:
            class_std_flat = None

        for i in range(class_difference):
            if method in ("average", "noisedaverage"):
                sample = data_class.sample(n=3, replace=True)
                name_up = f"up_{k}_{'_'.join(sample.index)}"
                k += 1
                profile_up = sample.median(axis=0).to_frame().transpose()

                if method == "noisedaverage":
                    nv = np.random.normal(
                        profile_up.values.flatten(),
                        noise_stds * class_std_flat,
                        size=profile_up.shape,
                    )
                    nv = np.where(nv > 1, 1, np.where(nv < 0, 0, nv))
                    profile_up = pd.DataFrame(nv, columns=profile_up.columns)
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
    # TODO: seems unnecessary?!
    fract_marker_up = fract_marker_up.sample(frac=1)

    assert fract_marker_up["class"].value_counts().nunique() == 1

    return fract_marker_up


def MOP_exec(
    fract_full: dict[str, pd.DataFrame],
    fract_marker: dict[str, pd.DataFrame],
    fract_test: dict[str, pd.DataFrame],
    nn_params: NeuralNetworkParametersModel,
    max_processes: int = 1,
) -> dict[str, ConditionPredictionModel]:
    """Perform multi-organelle prediction while showing a progress dialog.

    :param fract_full: dictionary of full profiles
    """
    import FreeSimpleGUI as sg

    def update_progress(
        progress_bars: dict[str, sg.ProgressBar], progress_queue: mp.Queue
    ):
        """Update the progress bars based on the progress queue."""
        while True:
            try:
                condition_id, round_id, percent, task = (
                    progress_queue.get_nowait()
                )
                finished[condition_id][round_id] = percent
                progress_bars[condition_id].update_bar(
                    int(sum(finished[condition_id].values()))
                )
                status_labels[condition_id].update(task)
            except queue.Empty:
                break

    # condition_id -> progress bar
    progress_bars = {}
    status_labels = {}
    # condition_id -> round_id -> percentage
    finished = {condition_id: {} for condition_id in fract_full}
    # maximum length of condition_ids for text box alignment
    text_maxlen = max(map(len, fract_full))
    layout = []
    for condition_id in fract_full:
        progress_bar = sg.ProgressBar(
            max_value=nn_params.rounds * 100,
            orientation="h",
            size=(20, 10),
            key=condition_id,
            expand_x=True,
            expand_y=True,
        )
        layout.append([sg.Text(condition_id, s=text_maxlen), progress_bar])
        progress_bars[condition_id] = progress_bar
        status_labels[condition_id] = sg.Text("")
        layout.append([status_labels[condition_id]])

    window = sg.Window(
        "Progress", layout, finalize=True, resizable=True, modal=True
    )
    window.set_cursor("watch")

    # run actual prediction in a separate thread to allow for progress updates
    manager = mp.Manager()
    progress_queue = manager.Queue()
    result_queue = mp.Queue()
    worker_thread = threading.Thread(
        target=multi_organelle_prediction,
        args=(
            fract_full,
            fract_marker,
            fract_test,
            nn_params,
            max_processes,
            progress_queue,
            result_queue,
        ),
        daemon=True,
    )
    worker_thread.start()

    while True:
        window.read(timeout=100)
        update_progress(
            progress_bars,
            progress_queue,
        )
        try:
            result = result_queue.get_nowait()
            break
        except queue.Empty:
            pass

    worker_thread.join()
    window.close()

    return result


def multi_organelle_prediction(
    fract_full: dict[str, pd.DataFrame],
    fract_marker: dict[str, pd.DataFrame],
    fract_test: dict[str, pd.DataFrame],
    nn_params: NeuralNetworkParametersModel,
    max_processes: int = 1,
    progress_queue: mp.Queue = None,
    result_queue: mp.Queue = None,
) -> dict[str, ConditionPredictionModel]:
    """Perform multi-organelle prediction.

    :param fract_full: dictionary of full profiles
    :param max_processes: The maximum number of processes to use.
    :param progress_queue: A queue to report progress
        `tuple[condition_id: str, round_id: str, percent_done: int|float,
        task: str]`.
    :param result_queue: A queue to report the result.
    """

    # prepare data structures for each round / condition
    conditions = list(fract_full.keys())
    class_predictions = {
        condition: ConditionPredictionModel(
            condition_id=condition,
            classes=fract_marker[condition]["class"].unique().tolist(),
        )
        for condition in conditions
    }

    # parallel execution
    args_list = [
        (
            condition,
            i_round,
            class_predictions[condition],
            fract_full[condition],
            fract_marker[condition],
            fract_test[condition],
            nn_params,
            logger,
            progress_queue,
        )
        for condition in conditions
        for i_round in range(1, nn_params.rounds + 1)
    ]

    def on_task_done(result):
        condition, round_id, round_data = result
        class_predictions[condition].round_results[round_id] = round_data

    if max_processes > 1:
        ctx = get_mp_ctx()
        with ctx.Pool(processes=max_processes) as pool:
            results = [
                pool.apply_async(
                    execute_round_wrapper, (args,), callback=on_task_done
                )
                for args in args_list
            ]

            for result in results:
                result.get()
    else:
        for args in args_list:
            result = execute_round_wrapper(args)
            on_task_done(result)

    if result_queue:
        result_queue.put(class_predictions)

    return class_predictions


def execute_round_wrapper(args):
    """Parallelization wrapper for `execute_round`."""
    (
        condition,
        i_round,
        class_predictions,
        fract_full,
        fract_marker,
        fract_test,
        nn_params,
        logger,
        progress_queue,
    ) = args
    logger.info(
        f"Executing round {i_round}/{nn_params.rounds} "
        f"for condition {condition}..."
    )
    round_id = f"ROUND_{i_round}"
    log_prefix = f"[{condition} {round_id}]"
    sub_logger = logger.getChild(log_prefix)
    sub_logger.addFilter(PrefixFilter(log_prefix))
    round_data = execute_round(
        class_predictions,
        fract_full,
        fract_marker,
        fract_test,
        nn_params,
        sub_logger,
        round_id,
        keras_proj_id=f"Classifier_{condition}_{i_round}",
        progress_queue=progress_queue,
    )
    return condition, round_id, round_data


def execute_round(
    class_predictions: ConditionPredictionModel,
    fract_full: pd.DataFrame,
    fract_marker: pd.DataFrame,
    fract_test: pd.DataFrame,
    nn_params: NeuralNetworkParametersModel,
    logger: logging.Logger,
    round_id: str,
    keras_proj_id: str,
    progress_queue: mp.Queue = None,
) -> TrainingRoundModel:
    """Perform a single round of training and prediction.

    :param fract_marker: The marker profiles.
    :param fract_test: The profiles of species with unknown classes.
    """

    result = TrainingRoundModel()

    if progress_queue:
        progress_queue.put(
            (class_predictions.condition_id, round_id, 0, "Upsampling...")
        )

    # upsample fractionation data
    if nn_params.upsampling:
        logger.info("Upsampling")
        fract_marker_up = upsample_condition(
            fract_marker,
            method=nn_params.upsampling_method,
            noise_stds=nn_params.upsampling_noise,
        )
        logger.info("upsampling done!")
    else:
        fract_marker_up = copy.deepcopy(fract_marker)

    if progress_queue:
        progress_queue.put(
            (class_predictions.condition_id, round_id, 1, "SVM prediction...")
        )

    logger.info("Performing single prediction ...")
    _, svm_marker, _, result.w_full_prob_df = single_prediction(
        x_train_up_df=fract_marker_up.drop(columns=["class"]),
        y_train_up_df=fract_marker_up["class"],
        fract_full=fract_full,
        fract_marker=fract_marker,
        fract_test=fract_test,
    )

    if nn_params.svm_filter:
        logger.info("Applying SVM filter...")
        # Remove the markers that are not predicted correctly by the SVM
        #  and upsample the rest
        fract_marker_filtered = fract_marker[
            svm_marker["class"] == svm_marker["svm_prediction"]
        ]
        logger.info("Upsampling after SVM-filtering...")
        fract_marker_up = upsample_condition(
            fract_marker_filtered,
            method=nn_params.upsampling_method,
            noise_stds=nn_params.upsampling_noise,
        )
        logger.info("SVM filtering done.")

    fract_unmixed_up = pd.concat(
        [
            fract_marker_up.drop("class", axis=1),
            # One-hot encode the classes
            pd.get_dummies(fract_marker_up["class"]),
        ],
        axis=1,
    )

    if progress_queue:
        progress_queue.put(
            (class_predictions.condition_id, round_id, 5, "Mixing profiles...")
        )

    if nn_params.mixed_part == "none":
        fract_mixed_up = copy.deepcopy(fract_unmixed_up)
    else:
        logger.info("Mixing profiles...")
        mix_steps = [
            i / nn_params.mixed_part for i in range(1, nn_params.mixed_part)
        ]
        fract_mixed_up = mix_profiles(
            fract_marker_up,
            fract_unmixed_up,
            mix_steps,
            nn_params.mixed_batch,
        )
        logger.info("mixing done!")

    if progress_queue:
        progress_queue.put(
            (class_predictions.condition_id, round_id, 10, "Training...")
        )

    def status_callback(percent_done: float, task: str):
        if progress_queue:
            # renormalize
            percent_done = percent_done / 100 * 90 + 10
            progress_queue.put(
                (class_predictions.condition_id, round_id, percent_done, task)
            )

    with stdout_to_logger(logger, logging.DEBUG):
        result.z_full_df = multi_predictions(
            x_full_df=fract_full.drop(columns=["class"]),
            x_train_df=fract_mixed_up.drop(columns=class_predictions.classes),
            y_train_df=fract_mixed_up[class_predictions.classes],
            nn_params=nn_params,
            logger=logger,
            keras_proj_id=keras_proj_id,
            status_callback=status_callback,
        )

    if progress_queue:
        progress_queue.put(
            (class_predictions.condition_id, round_id, 100, "done")
        )

    return result


def multi_predictions(
    x_full_df: pd.DataFrame,
    x_train_df: pd.DataFrame,
    y_train_df: pd.DataFrame,
    nn_params: NeuralNetworkParametersModel,
    logger: logging.Logger,
    keras_proj_id: str,
    status_callback: Callable[[float, str], None] = lambda *args,
    **kwargs: None,
) -> pd.DataFrame:
    """Perform multi organelle predictions.

    1. Perform hyperparameter tuning for the neural network.
    2. Re-train the neural network with the best hyperparameters and predict
    the full dataset.
    3. Return the averaged predictions

    :param x_full_df: The profiles for which to predict the classes
        (species × fraction).
    :param x_train_df: The training profiles (species × fraction).
    :param y_train_df: The classes for the training profiles (species × class).
    :param nn_params: The neural network parameters.
    :param logger: A logger.
    :param keras_proj_id: The project ID for the keras tuner.
    :param status_callback: A callback to report the status.
    :return: The predicted classes for the full profiles. The average of the
        neural network predictions from all re-training/prediction rounds.
    """
    logger.info("Training classifier...")

    import keras_tuner as kt
    import tensorflow as tf

    from .classification_model import FNN_Classifier

    x_train = x_train_df.to_numpy(dtype=float)
    y_train = y_train_df.to_numpy(dtype=float)
    num_compartments = y_train.shape[1]
    num_fractions = x_train.shape[1]
    set_shapes = [num_fractions, num_compartments]

    # Set up hyperparameter tuning
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
        project_name=f"{time}_{keras_proj_id}",
    )
    # Tune the hyperparameters
    status_callback(0, "Hyperparameter tuning...")
    # noinspection PyUnresolvedReferences
    stop_early = tf.keras.callbacks.EarlyStopping(
        monitor="val_loss", patience=5
    )
    tuner.search(
        x_train,
        y_train,
        epochs=nn_params.NN_epochs,
        validation_split=0.2,
        callbacks=[stop_early],
    )
    logger.info("Hyperparameter tuning done!")

    best_model = tuner.get_best_models(num_models=1)[0]
    best_hp = tuner.get_best_hyperparameters(num_trials=1)[0].values

    # Re-train the model with the best hyperparameters and predict the full
    #  dataset.
    #  (For the first subround, we use the results of the hyperparameter tuning
    #  directly).
    z_full_arrays = [best_model.predict(x_full_df.values)]

    for i_subround in range(2, nn_params.subrounds + 1):
        status_callback(
            50 + i_subround / nn_params.subrounds * 50,
            f"Training round {i_subround}/{nn_params.subrounds}...",
        )
        logger.info(
            f"Training classifier for {i_subround}/{nn_params.subrounds}..."
        )

        fixed_model = FNN_Classifier(
            fixed_hp=best_hp, set_shapes=set_shapes, nn_params=nn_params
        ).build(None)
        fixed_model.fit(
            x_train,
            y_train,
            epochs=nn_params.NN_epochs,
            validation_split=0.2,
            callbacks=[stop_early],
        )

        z_full_arrays.append(fixed_model.predict(x_full_df.values))

    # free memory
    tf.keras.backend.clear_session()
    gc.collect()

    # average predictions
    z_full_df = pd.DataFrame(
        np.stack(z_full_arrays).mean(axis=0),
        index=x_full_df.index,
        columns=y_train_df.columns,
    )
    return z_full_df


def mix_profiles(
    fract_marker_up: pd.DataFrame,
    fract_unmixed_up: pd.DataFrame,
    mix_steps: list[float],
    mixed_batch: float,
) -> pd.DataFrame:
    """Create mixed profiles.

    Complement `fract_unmixed_up` with artificially mixed marker profiles.

    For each combination of classes/compartments, create mixed profiles,
    by pairwise combination of different marker profiles (`fract_marker_up`)
    according to the given fractions in `mix_steps`. Then subsample the mixed
    profiles.
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
        profiles_own = fract_marker_up.loc[
            fract_marker_up["class"] == comb[0]
        ].drop(columns=["class"])
        profiles_other = fract_marker_up.loc[
            fract_marker_up["class"] == comb[1]
        ].drop(columns=["class"])

        new_index = [
            f"{i}_{j}"
            for i, j in zip(profiles_own.index, profiles_other.index)
        ]
        for part in mix_steps:
            new_index_part = [
                f"{i + cur}_{value}" for i, value in enumerate(new_index)
            ]
            # Actual mixing
            own_part = profiles_own.multiply(part)
            other_part = profiles_other.multiply(1 - part)

            own_part.index = new_index_part
            other_part.index = new_index_part

            profiles_mixed = own_part + other_part

            # Add columns of class ratios
            for classname in class_list:
                if classname == comb[0]:
                    profiles_mixed[classname] = part
                elif classname == comb[1]:
                    profiles_mixed[classname] = 1 - part
                else:
                    profiles_mixed[classname] = 0.0

            # subsample mixed profiles
            profiles_mixed = profiles_mixed.sample(frac=mixed_batch)

            fract_mixed_up = pd.concat(
                [fract_mixed_up, profiles_mixed]
            ).sample(frac=1)
            cur += len(profiles_mixed)
    return fract_mixed_up


def single_prediction(
    x_train_up_df: pd.DataFrame,
    y_train_up_df: pd.Series,
    fract_full: pd.DataFrame,
    fract_marker: pd.DataFrame,
    fract_test: pd.DataFrame,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Perform single-class (single-compartment) prediction.

    Train Support Vector Machine (SVM) classifier and predict the classes.

    :param x_train_up_df: The upsampled training profiles.
    :param y_train_up_df: The classes for the training profiles.
    :param fract_full: Profiles with known and unknown classes.
    :param fract_marker: The marker profiles.
    :param fract_test: The profiles of species with unknown classes.
    """
    # train classifier on the upsampled data
    clf = svm.SVC(kernel="rbf", probability=True)
    clf.fit(x_train_up_df, y_train_up_df)

    # predict for the full dataset
    x_full = fract_full.drop(columns=["class"])
    w_full = clf.predict(x_full).tolist()
    w_full_prob = clf.predict_proba(x_full).max(axis=1)
    w_full_prob_df = pd.DataFrame(
        {
            "SVM_winner": w_full,
            "SVM_prob": w_full_prob,
        },
        index=x_full.index,
    )

    # evaluate on the (non-artificial) marker profiles
    x_train = fract_marker.drop(columns=["class"])
    W_train = fract_marker["class"]
    w_train = clf.predict(x_train).tolist()
    w_train_prob = clf.predict_proba(x_train).max(axis=1)
    svm_marker = pd.DataFrame(
        {
            "class": W_train,
            "svm_prediction": w_train,
            "svm_probability": w_train_prob,
        },
        index=x_train.index,
    )

    # compute stats on the marker classification
    # TODO: unused
    confusion = pd.DataFrame(
        confusion_matrix(W_train, w_train, labels=list(clf.classes_)),
        index=clf.classes_,
        columns=clf.classes_,
    )
    accuracy = accuracy_score(W_train, w_train)
    precision = precision_score(W_train, w_train, average="macro")
    recall = recall_score(W_train, w_train, average="macro")
    f1 = f1_score(W_train, w_train, average="macro")
    svm_metrics = {
        "confusion": confusion,
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }

    # predict species with unknown localization
    # TODO: unused
    x_test = fract_test.drop(columns=["class"])
    w_test = clf.predict(x_test).tolist()
    w_test_prob = clf.predict_proba(x_test).max(axis=1)
    svm_test = pd.DataFrame(
        {
            "svm_prediction": w_test,
            "svm_probability": w_test_prob,
        },
        index=x_test.index,
    )

    return svm_metrics, svm_marker, svm_test, w_full_prob_df
