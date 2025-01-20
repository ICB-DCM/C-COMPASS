"""Single organelle prediction."""
# FIXME: this module is currently unused

import copy

import matplotlib.pyplot as plt
import pandas as pd
from sklearn import svm
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
)

clf = svm.SVC(kernel="rbf", probability=True)


def get_cmap(n, name="hsv"):
    return plt.cm.get_cmap(name, n)


def plot_compressed(data, classes, name):
    cmap = get_cmap(len(classes))
    fig, ax = plt.subplots()
    for i in range(len(classes)):
        organelle = classes[i]
        data_class = data.loc[data["class"] == organelle]
        current_color = cmap(i)
        if organelle == "none":
            current_color = "grey"
        ax.scatter(
            data_class[0],
            data_class[1],
            s=0.2,
            color=current_color,
            label=organelle,
        )
    ax.legend(fontsize="x-small")
    plt.title(name)
    plt.show()


def create_wxyz(
    fract_conditions,
    fract_marker,
    fract_marker_up,
    fract_test,
    fract_full,
    fract_full_up,
    fract_mixed_up,
):
    # split_size = 0.8

    conditions = [x for x in fract_conditions if x != "[KEEP]"]
    learning_wxyz = {}

    # df_randomized = df.sample(frac=1)

    for condition in conditions:
        fract_marker[condition] = fract_marker[condition].sample(frac=1)
        fract_marker_up[condition] = fract_marker_up[condition].sample(frac=1)
        fract_mixed_up[condition] = fract_mixed_up[condition].sample(frac=1)

        # Calculate the split index
        # split_index = int(split_size * len(fract_marker_up[condition]))

        # Split the DataFrame
        # fract_train_up_con = fract_marker_up[condition][:split_index]
        # fract_val_up_con = fract_marker_up[condition][split_index:]

        classes = list(set(fract_marker[condition]["class"]))
        learning_wxyz[condition] = {}

        W_train = list(fract_marker[condition]["class"])
        learning_wxyz[condition]["W_train"] = W_train

        W_train_up = list(fract_marker_up[condition]["class"])
        # W_train_up = list(fract_train_up_con['class'])
        learning_wxyz[condition]["W_train_up"] = W_train_up

        # W_val_up = list(fract_val_up_con['class'])
        # learning_wxyz[condition]['W_val_up'] = W_val_up

        W_full = list(fract_full[condition]["class"])
        learning_wxyz[condition]["W_full"] = W_full

        W_full_up = list(fract_full_up[condition]["class"])
        learning_wxyz[condition]["W_full_up"] = W_full_up

        x_train = (
            fract_marker[condition]
            .drop(columns=["class"])
            .to_numpy(dtype=float)
        )
        learning_wxyz[condition]["x_train"] = x_train

        x_train_up = (
            fract_marker_up[condition]
            .drop(columns=["class"])
            .to_numpy(dtype=float)
        )
        # x_train_up = fract_train_up_con.drop(columns = ['class']).to_numpy(dtype = float)
        learning_wxyz[condition]["x_train_up"] = x_train_up

        # x_val_up = fract_val_up_con.drop(columns = ['class']).to_numpy(dtype = float)
        # learning_wxyz[condition]['x_val_up'] = x_val_up

        x_test = (
            fract_test[condition].drop(columns=["class"]).to_numpy(dtype=float)
        )
        learning_wxyz[condition]["x_test"] = x_test

        x_full = (
            fract_full[condition].drop(columns=["class"]).to_numpy(dtype=float)
        )
        learning_wxyz[condition]["x_full"] = x_full

        x_full_up = (
            fract_full_up[condition]
            .drop(columns=["class"])
            .to_numpy(dtype=float)
        )
        learning_wxyz[condition]["x_full_up"] = x_full_up

        x_train_mixed_up = (
            fract_mixed_up[condition]
            .drop(columns=classes)
            .to_numpy(dtype=float)
        )
        learning_wxyz[condition]["x_train_mixed_up"] = x_train_mixed_up

        Z_train = pd.get_dummies(fract_marker[condition]["class"]).to_numpy(
            dtype=float
        )
        learning_wxyz[condition]["Z_train"] = Z_train

        Z_train_up = pd.get_dummies(
            fract_marker_up[condition]["class"]
        ).to_numpy(dtype=float)
        # Z_train_up = pd.get_dummies(fract_train_up_con['class']).to_numpy(dtype = float)
        learning_wxyz[condition]["Z_train_up"] = Z_train_up

        # Z_val_up = pd.get_dummies(fract_val_up_con['class']).to_numpy(dtype = float)
        # learning_wxyz[condition]['Z_val_up'] = Z_val_up

        Z_train_mixed_up = fract_mixed_up[condition][classes].to_numpy(
            dtype=float
        )
        learning_wxyz[condition]["Z_train_mixed_up"] = Z_train_mixed_up

    return learning_wxyz


# def create_xZ (fract_full, fract_full_up, fract_marker, fract_marker_up, fract_test, fract_conditions):
#     conditions = [x for x in fract_conditions if x != '[KEEP]']
#     learning_xZ = {}
#     for condition in conditions:
#         learning_xZ[condition] = {}

#         x_train = fract_marker[condition].drop(columns = ['class']).to_numpy(dtype = float)
#         learning_xZ[condition]['x_train'] = x_train

#         x_train_up = fract_marker_up[condition].drop(columns = ['class']).to_numpy(dtype = float)
#         learning_xZ[condition]['x_train_up'] = x_train_up

#         x_test = fract_test[condition].drop(columns = ['class']).to_numpy(dtype = float)
#         learning_xZ[condition]['x_test'] = x_test

#         x_full = fract_full[condition].drop(columns = ['class']).to_numpy(dtype = float)
#         learning_xZ[condition]['x_full'] = x_full

#         x_full_up = fract_full_up[condition].drop(columns = ['class']).to_numpy(dtype = float)
#         learning_xZ[condition]['x_full_up'] = x_full_up

#         Z_train = list(fract_marker[condition]['class'])
#         learning_xZ[condition]['Z_train'] = Z_train

#         Z_train_up = list(fract_marker_up[condition]['class'])
#         learning_xZ[condition]['Z_train_up'] = Z_train_up

#         Z_full = list(fract_full[condition]['class'])
#         learning_xZ[condition]['Z_full'] = Z_full

#         Z_full_up = list(fract_full_up[condition]['class'])
#         learning_xZ[condition]['Z_full_up'] = Z_full_up
#     return learning_xZ


# def SOP_exec (fract_full, fract_full_up, fract_marker, fract_marker_up, fract_test, fract_conditions):


def SOP_exec(learning_wxyz, fract_conditions, fract_marker, fract_test):
    conditions = [x for x in fract_conditions if x != "[KEEP]"]
    svm_metrics = {}

    svm_marker = {}
    svm_test = {}

    for condition in conditions:
        svm_x_train = learning_wxyz[condition]["x_train"]
        svm_x_train_up = learning_wxyz[condition]["x_train_up"]
        svm_x_test = learning_wxyz[condition]["x_test"]
        svm_Y_train = learning_wxyz[condition]["W_train"]
        svm_Y_train_up = learning_wxyz[condition]["W_train_up"]

        clf.fit(svm_x_train_up, svm_Y_train_up)

        svm_y_train = list(clf.predict(svm_x_train))
        svm_y_test = list(clf.predict(svm_x_test))

        svm_y_train_prob = list(
            map(lambda x: max(x), list(clf.predict_proba(svm_x_train)))
        )
        svm_y_test_prob = list(
            map(lambda x: max(x), list(clf.predict_proba(svm_x_test)))
        )

        svm_confusion = pd.DataFrame(
            confusion_matrix(
                svm_Y_train, svm_y_train, labels=list(clf.classes_)
            ),
            index=clf.classes_,
            columns=clf.classes_,
        )
        svm_accuracy = accuracy_score(svm_Y_train, svm_y_train)
        svm_precision = precision_score(
            svm_Y_train, svm_y_train, average="macro"
        )
        svm_recall = recall_score(svm_Y_train, svm_y_train, average="macro")
        svm_f1 = f1_score(svm_Y_train, svm_y_train, average="macro")

        svm_marker[condition] = copy.deepcopy(fract_marker[condition])
        svm_marker[condition]["svm_prediction"] = svm_y_train
        svm_marker[condition]["svm_probability"] = svm_y_train_prob

        svm_test[condition] = copy.deepcopy(fract_test[condition])
        svm_test[condition]["svm_prediction"] = svm_y_test
        svm_test[condition]["svm_probability"] = svm_y_test_prob

        svm_metrics[condition] = {
            "confusion": svm_confusion,
            "accuracy": svm_accuracy,
            "precision": svm_precision,
            "recall": svm_recall,
            "f1": svm_f1,
        }

    return svm_marker, svm_test, svm_metrics


def plot_umaps(
    learning_wxyz,
    fract_marker,
    fract_marker_up,
    fract_test,
    fract_full,
    fract_full_up,
):
    import umap.umap_ as umap

    umaps = {}

    for condition in learning_wxyz:
        umaps[condition] = {}

        classes = list(set(learning_wxyz[condition]["W_train"]))

        x_train = learning_wxyz[condition]["x_train"]
        x_train_up = learning_wxyz[condition]["x_train_up"]
        x_test = learning_wxyz[condition]["x_test"]
        x_full = learning_wxyz[condition]["x_full"]
        x_full_up = learning_wxyz[condition]["x_full_up"]

        classes_dict = {}
        i = 0
        for cl in classes:
            classes_dict[cl] = i
            i = i + 1
        classes_dict["none"] = i

        classes.append("none")

        reducer = umap.UMAP()

        train_umap = reducer.fit_transform(x_train)
        train_umap_df = pd.DataFrame(train_umap)
        train_umap_df["class"] = list(fract_marker[condition]["class"])
        plot_compressed(train_umap_df, classes, "UMAP_train_" + condition)

        train_up_umap = reducer.fit_transform(x_train_up)
        train_up_umap_df = pd.DataFrame(train_up_umap)
        train_up_umap_df["class"] = list(fract_marker_up[condition]["class"])
        plot_compressed(
            train_up_umap_df, classes, "UMAP_train_up_" + condition
        )

        test_umap = reducer.fit_transform(x_test)
        test_umap_df = pd.DataFrame(test_umap)
        test_umap_df["class"] = list(
            fract_test[condition].fillna("none")["class"]
        )
        plot_compressed(test_umap_df, classes, "UMAP_test_" + condition)

        full_umap = reducer.fit_transform(x_full)
        full_umap_df = pd.DataFrame(full_umap)
        full_umap_df["class"] = list(
            fract_full[condition].fillna("none")["class"]
        )
        plot_compressed(full_umap_df, classes, "UMAP_full_" + condition)

        full_up_umap = reducer.fit_transform(x_full_up)
        full_up_umap_df = pd.DataFrame(full_up_umap)
        full_up_umap_df["class"] = list(
            fract_full_up[condition].fillna("none")["class"]
        )
        plot_compressed(full_up_umap_df, classes, "UMAP_full_up_" + condition)

        umaps[condition]["svm_train_umap"] = train_umap_df
        umaps[condition]["svm_train_up_umap"] = train_up_umap_df
        umaps[condition]["svm_test_umap"] = test_umap_df
        umaps[condition]["svm_full_umap"] = full_umap_df
        umaps[condition]["svm_full_up_umap"] = full_up_umap_df

    return
