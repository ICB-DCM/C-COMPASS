"""Multi-organelle prediction statistics."""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import f, ttest_ind


def reduce_to_single_column(row):
    if row.nunique() == 1:
        return row.iloc[0]
    else:
        return np.nan


def calculate_mean(row):
    if pd.isnull(row["SVM_winner"]):
        return np.nan
    row_values = row.dropna()[
        1:
    ]  # Exclude the first column which is 'reduced_column'
    return np.mean(row_values)


def most_frequent_or_nan(row):
    counts = row.value_counts()
    # If the row is empty, return np.nan
    if counts.empty:
        return np.nan
    # If there's only one unique value in the row, return that value
    elif len(counts) == 1:
        return counts.idxmax()
    # If the two most frequent values occur the same number of times, return np.nan
    elif counts.iloc[0] == counts.iloc[1]:
        return np.nan
    else:
        return counts.idxmax()


def safe_sum(lst):
    return np.nansum(lst)


# def most_frequent_or_nan(row):
#     counts = row.value_counts()
#     # Check if the row is empty or if the Series has less than 2 elements
#     if len(counts) < 2:
#         return np.nan
#     # Check if the two most frequent values occur the same number of times
#     elif counts.iloc[0] == counts.iloc[1]:
#         return np.nan
#     else:
#         return counts.idxmax()

# def most_frequent_or_nan(row):
#     counts = row.value_counts()
#     if len(counts) == 0 or counts.iloc[0] == counts.iloc[1]:
#         return np.nan
#     return counts.idxmax()


def stats_exec3(
    mode,
    learning_xyz,
    NN_params,
    fract_data,
    tp_data,
    fract_marker,
    fract_conditions,
    fract_std,
):
    conditions = [x for x in fract_conditions if x != "[KEEP]"]
    results = {}

    for condition in conditions:
        subcons = [x for x in learning_xyz if x.startswith(condition + "_")]

        if mode == "deep":
            combined_index = tp_data[condition].index
        elif mode == "rough":
            empty_index = []
            df = pd.DataFrame(index=empty_index)
            combined_index = df.index
        for subcon in subcons:
            combined_index = combined_index.union(
                fract_data["class"][subcon].index
            )

        # combined_index = tp_data[condition].index
        # for subcon in subcons:
        #     combined_index = combined_index.union(fract_data['class'][subcon].index)

        results[condition] = {}
        results[condition]["metrics"] = pd.DataFrame(index=combined_index)

        TPA_list = []

        ## add TPA:
        if mode == "deep":
            TPA_list = []
            tp_nontrans = tp_data[condition].map(lambda x: 2**x)
            for replicate in tp_data[condition]:
                TPA_list.append(tp_nontrans[replicate])
            combined_TPA = pd.concat(TPA_list, axis=1)
            results[condition]["metrics"]["TPA"] = combined_TPA.mean(axis=1)
            results[condition]["metrics"] = results[condition]["metrics"].loc[
                ~results[condition]["metrics"].index.duplicated(keep="first")
            ]
        elif mode == "rough":
            pass

        ## add marker:
        results[condition]["metrics"]["marker"] = np.nan
        # df_2 = df_2[~df_2.index.duplicated(keep='first')]

        for subcon in subcons:
            marker_df = learning_xyz[subcon]["W_train_df"][
                ~learning_xyz[subcon]["W_train_df"].index.duplicated(
                    keep="first"
                )
            ]
            results[condition]["metrics"]["marker"] = results[condition][
                "metrics"
            ]["marker"].fillna(marker_df)

        ## add SVM results:
        results[condition]["SVM"] = {}
        results[condition]["SVM"]["winner_combined"] = pd.DataFrame(
            index=results[condition]["metrics"].index
        )
        results[condition]["SVM"]["prob_combined"] = pd.DataFrame(
            index=results[condition]["metrics"].index
        )
        for subcon in subcons:
            print(subcon)
            learning_xyz[subcon]["w_full_combined"] = pd.DataFrame(
                index=learning_xyz[subcon]["x_full_df"].index
            )
            learning_xyz[subcon]["w_full_prob_combined"] = pd.DataFrame(
                index=learning_xyz[subcon]["x_full_df"].index
            )

            for roundn in learning_xyz[subcon]["w_full_prob_df"]:
                learning_xyz[subcon]["w_full_combined"] = pd.merge(
                    learning_xyz[subcon]["w_full_combined"],
                    learning_xyz[subcon]["w_full_prob_df"][roundn][
                        "SVM_winner"
                    ],
                    left_index=True,
                    right_index=True,
                    how="left",
                )
                learning_xyz[subcon]["w_full_combined"] = learning_xyz[subcon][
                    "w_full_combined"
                ].loc[
                    ~learning_xyz[subcon]["w_full_combined"].index.duplicated(
                        keep="first"
                    )
                ]
                learning_xyz[subcon]["w_full_combined"].rename(
                    columns={"SVM_winner": roundn + "_SVM_winner"},
                    inplace=True,
                )
                learning_xyz[subcon]["w_full_prob_combined"] = pd.merge(
                    learning_xyz[subcon]["w_full_prob_combined"],
                    learning_xyz[subcon]["w_full_prob_df"][roundn]["SVM_prob"],
                    left_index=True,
                    right_index=True,
                    how="left",
                )
                learning_xyz[subcon]["w_full_prob_combined"] = learning_xyz[
                    subcon
                ]["w_full_prob_combined"].loc[
                    ~learning_xyz[subcon][
                        "w_full_prob_combined"
                    ].index.duplicated(keep="first")
                ]
                learning_xyz[subcon]["w_full_prob_combined"].rename(
                    columns={"SVM_prob": roundn + "_SVM_prob"}, inplace=True
                )

            SVM_equal = learning_xyz[subcon]["w_full_combined"].apply(
                lambda row: row.nunique() == 1, axis=1
            )
            learning_xyz[subcon]["w_full_combined"]["SVM_winner"] = np.where(
                SVM_equal,
                learning_xyz[subcon]["w_full_combined"].iloc[:, 0],
                np.nan,
            )
            learning_xyz[subcon]["w_full_prob_combined"]["SVM_prob"] = (
                learning_xyz[subcon]["w_full_prob_combined"].mean(axis=1)
            )

            learning_xyz[subcon]["w_full_combined"] = pd.merge(
                learning_xyz[subcon]["w_full_combined"],
                learning_xyz[subcon]["w_full_prob_combined"][["SVM_prob"]],
                left_index=True,
                right_index=True,
                how="left",
            )
            learning_xyz[subcon]["w_full_combined"] = learning_xyz[subcon][
                "w_full_combined"
            ].loc[
                ~learning_xyz[subcon]["w_full_combined"].index.duplicated(
                    keep="first"
                )
            ]
            learning_xyz[subcon]["w_full_combined"].loc[
                learning_xyz[subcon]["w_full_combined"]["SVM_winner"].isna(),
                "SVM_prob",
            ] = np.nan

            results[condition]["SVM"]["winner_combined"] = pd.merge(
                results[condition]["SVM"]["winner_combined"],
                learning_xyz[subcon]["w_full_combined"]["SVM_winner"],
                left_index=True,
                right_index=True,
                how="left",
            )
            results[condition]["SVM"]["winner_combined"] = results[condition][
                "SVM"
            ]["winner_combined"].loc[
                ~results[condition]["SVM"]["winner_combined"].index.duplicated(
                    keep="first"
                )
            ]
            results[condition]["SVM"]["prob_combined"] = pd.merge(
                results[condition]["SVM"]["prob_combined"],
                learning_xyz[subcon]["w_full_combined"]["SVM_prob"],
                left_index=True,
                right_index=True,
                how="left",
            )
            results[condition]["SVM"]["prob_combined"] = results[condition][
                "SVM"
            ]["prob_combined"].loc[
                ~results[condition]["SVM"]["prob_combined"].index.duplicated(
                    keep="first"
                )
            ]

        SVM_equal = results[condition]["SVM"]["winner_combined"].apply(
            lambda row: row.nunique() == 1, axis=1
        )
        # SVM_equal =
        SVM_major = results[condition]["SVM"]["winner_combined"].apply(
            most_frequent_or_nan, axis=1
        )
        SVM_major.name = "SVM_subwinner"
        results[condition]["metrics"]["SVM_winner"] = np.where(
            SVM_equal,
            results[condition]["SVM"]["winner_combined"].iloc[:, 0],
            np.nan,
        )
        results[condition]["metrics"] = pd.merge(
            results[condition]["metrics"],
            SVM_major,
            left_index=True,
            right_index=True,
            how="left",
        )
        prob_means = results[condition]["SVM"]["prob_combined"].mean(axis=1)
        results[condition]["metrics"]["SVM_prob"] = np.nan
        results[condition]["metrics"].loc[
            results[condition]["metrics"]["SVM_winner"].notna(), "SVM_prob"
        ] = prob_means

        ## add CClist:
        for subcon in subcons:
            stacked_arrays = np.stack(
                list(learning_xyz[subcon]["z_full"].values())
            )
            learning_xyz[subcon]["z_full_mean"] = np.mean(
                stacked_arrays, axis=0
            )
            learning_xyz[subcon]["z_full_mean_df"] = pd.DataFrame(
                learning_xyz[subcon]["z_full_mean"],
                index=learning_xyz[subcon]["x_full_df"].index,
                columns=learning_xyz[subcon]["Z_train_df"].columns,
            )

        con_frames = []
        for subcon in subcons:
            con_frames.append(learning_xyz[subcon]["z_full_mean_df"])

        classnames = []
        for subcon in subcons:
            for classname in learning_xyz[subcon]["classes"]:
                classnames.append(classname)
        classnames = list(set(classnames))

        for classname in classnames:
            CC_list = pd.DataFrame(index=combined_index)
            for subcon in subcons:
                CC_list = pd.merge(
                    CC_list,
                    learning_xyz[subcon]["z_full_mean_df"][classname],
                    left_index=True,
                    right_index=True,
                    how="left",
                )

            CC_list["CClist_" + classname] = CC_list.apply(
                lambda row: row.tolist(), axis=1
            )

            results[condition]["metrics"] = pd.merge(
                results[condition]["metrics"],
                CC_list["CClist_" + classname],
                left_index=True,
                right_index=True,
                how="left",
            )
            results[condition]["metrics"] = results[condition]["metrics"].loc[
                ~results[condition]["metrics"].index.duplicated(keep="first")
            ]

        ## add CC:
        for classname in classnames:
            # results[condition]['metrics']['CC_'+classname] = results[condition]['metrics']['CClist_'+classname].apply(lambda x: np.mean(x) if x else np.nan)
            results[condition]["metrics"]["CC_" + classname] = results[
                condition
            ]["metrics"]["CClist_" + classname].apply(
                lambda x: np.nanmean(x) if x else np.nan
            )
        cc_cols = [
            col
            for col in results[condition]["metrics"].columns
            if col.startswith("CC_")
        ]
        cc_sums = results[condition]["metrics"][cc_cols].sum(
            axis=1, skipna=True
        )
        # cc_sums = results[condition]['metrics'][cc_cols].map(safe_sum)
        results[condition]["metrics"][cc_cols] = results[condition]["metrics"][
            cc_cols
        ].div(cc_sums, axis=0)

        # cc_sums = results[condition]['metrics'][cc_columns].sum(axis=1, skipna=True)

        # results[condition]['metrics'][cc_columns] = results[condition]['metrics'][cc_columns].div(cc_sums, axis=0)
        # results[condition]['metrics'][cc_columns] = results[condition]['metrics'][cc_columns].div(cc_sums, axis=0)

        # pass

        ## add NN_winner:
        cc_columns = results[condition]["metrics"][
            [
                col
                for col in results[condition]["metrics"].columns
                if col.startswith("CC_")
            ]
        ]
        max_col = cc_columns.idxmax(axis=1)
        # print(max_col)
        # print(max_col.dtype)
        # print(max_col.head())
        max_col = max_col.astype(str)
        results[condition]["metrics"]["NN_winner"] = max_col.str.replace(
            "CC_", ""
        )

        ## add CA:
        if mode == "deep":
            results[condition]["metrics"]["CA_relevant"] = "no"
            results[condition]["class_abundance"] = {}
            for classname in classnames:
                # for classname in [classname for classname in results[condition]['metrics']['SVM_winner'] if not pd.isnull(classname)]:
                results_class = results[condition]["metrics"][
                    (results[condition]["metrics"]["NN_winner"] == classname)
                    & (~results[condition]["metrics"]["TPA"].isnull())
                ]
                results[condition]["metrics"].loc[
                    results_class.index, "CA_relevant"
                ] = "yes"
                results[condition]["class_abundance"][classname] = {}
                results[condition]["class_abundance"][classname]["CA"] = (
                    np.median(results_class["TPA"])
                )
                results[condition]["class_abundance"][classname]["count"] = (
                    len(results_class)
                )
        elif mode == "rough":
            pass

        ## add NN_winner:
        # columns_to_compare = [col for col in results[condition]['metrics'].columns if col.startswith('CC_')]
        # results[condition]['metrics']['NN_winner'] = results[condition]['metrics'][columns_to_compare].idxmax(axis=1)
        # results[condition]['metrics']['NN_winner'] = results[condition]['metrics']['NN_winner'].str.replace('^CC_', '', regex=True)

        # add fCC:
        for class_act in classnames:
            nonmarker_z = results[condition]["metrics"].loc[
                (results[condition]["metrics"]["marker"] != class_act)
                & (results[condition]["metrics"]["marker"].isna() == False)
            ][["CC_" + class_act]]
            thresh = np.percentile(
                nonmarker_z["CC_" + class_act].tolist(),
                NN_params["reliability"],
            )
            results[condition]["metrics"]["fCC_" + class_act] = results[
                condition
            ]["metrics"]["CC_" + class_act]
            results[condition]["metrics"].loc[
                results[condition]["metrics"]["fCC_" + class_act] < thresh,
                "fCC_" + class_act,
            ] = 0.0

        fcc_cols = [
            col
            for col in results[condition]["metrics"]
            if col.startswith("fCC_")
        ]
        fcc_sums = results[condition]["metrics"][fcc_cols].sum(axis=1)
        # fcc_sums[fcc_sums == 0] = 1
        results[condition]["metrics"][fcc_cols] = results[condition][
            "metrics"
        ][fcc_cols].div(fcc_sums, axis=0)

        ## add fNN_winner:
        fcc_columns = results[condition]["metrics"][
            [
                col
                for col in results[condition]["metrics"].columns
                if col.startswith("fCC_")
            ]
        ]
        fmax_col = fcc_columns.idxmax(axis=1)
        fmax_col = fmax_col.astype(str)
        results[condition]["metrics"]["fNN_winner"] = fmax_col.str.replace(
            "fCC_", ""
        )

        ## add nCClist:
        if mode == "deep":
            for classname in classnames:
                results[condition]["metrics"]["nCClist_" + classname] = (
                    results[condition]["metrics"]["CClist_" + classname].apply(
                        lambda lst: [
                            x
                            * results[condition]["class_abundance"][classname][
                                "CA"
                            ]
                            if not np.isnan(x)
                            else np.nan
                            for x in lst
                        ]
                    )
                )

            # nCClist_df = results[condition]['metrics'][]

            ## add nCC:
            for classname in classnames:
                results[condition]["metrics"]["nCC_" + classname] = (
                    results[condition]["metrics"]["fCC_" + classname]
                    * results[condition]["class_abundance"][classname]["CA"]
                )
            # normalize:
            nCC_cols = [
                col
                for col in results[condition]["metrics"]
                if col.startswith("nCC_")
            ]
            nCC_sums = results[condition]["metrics"][nCC_cols].sum(axis=1)
            nCC_sums[nCC_sums == 0] = 1
            results[condition]["metrics"][nCC_cols] = results[condition][
                "metrics"
            ][nCC_cols].div(nCC_sums, axis=0)

            ## add CPA
            for classname in classnames:
                results[condition]["metrics"]["CPA_" + classname] = (
                    results[condition]["metrics"]["CC_" + classname]
                    * results[condition]["metrics"]["TPA"]
                )
            for classname in classnames:
                results[condition]["metrics"]["nCPA_" + classname] = (
                    results[condition]["metrics"]["nCC_" + classname]
                    * results[condition]["metrics"]["TPA"]
                )
        elif mode == "rough":
            pass

        results[condition]["classnames"] = classnames

    return results


def stats_exec2(learning_xyz, NN_params, tp_data, fract_marker):
    conditions = [x for x in learning_xyz if x != "[KEEP]"]
    results = {}

    for condition in conditions:
        print(condition + "...")
        results[condition] = {}
        results[condition]["metrics"] = pd.DataFrame(
            index=learning_xyz[condition]["x_full_df"].index
        )

        # add TPA:
        TPA_list = []
        tp_nontrans = tp_data[condition].map(lambda x: 2**x)
        for replicate in tp_data[condition]:
            TPA_list.append(tp_nontrans[replicate])

        combined_TPA = pd.concat(TPA_list, axis=1)
        results[condition]["metrics"]["TPA"] = combined_TPA.mean(axis=1)

        # add SVM results:
        winner_list = []
        score_list = []
        for roundn in learning_xyz[condition]["w_full_prob_df"]:
            winner_list.append(
                learning_xyz[condition]["w_full_prob_df"][roundn]["SVM_winner"]
            )
            score_list.append(
                learning_xyz[condition]["w_full_prob_df"][roundn]["SVM_prob"]
            )
        combined_winner = pd.concat(winner_list, axis=1)
        results[condition]["metrics"]["SVM_winner"] = combined_winner.apply(
            reduce_to_single_column, axis=1
        )
        combined_score = pd.concat(score_list, axis=1)
        results[condition]["metrics"]["SVM_prob"] = combined_score.mean(axis=1)
        results[condition]["metrics"].loc[
            results[condition]["metrics"]["SVM_winner"].isna(), "SVM_prob"
        ] = np.nan

        # add CA:
        results[condition]["metrics"] = pd.merge(
            results[condition]["metrics"],
            fract_marker[condition]["class"],
            left_index=True,
            right_index=True,
            how="left",
        ).rename(columns={"class": "marker"})
        results[condition]["metrics"]["CA_relevant"] = "no"
        results[condition]["class_abundance"] = {}
        # classnames = []
        for classname in [
            classname
            for classname in results[condition]["metrics"]["SVM_winner"]
            if not pd.isnull(classname)
        ]:
            results_class = results[condition]["metrics"][
                (results[condition]["metrics"]["SVM_winner"] == classname)
                & (results[condition]["metrics"]["SVM_prob"] >= 0.90)
                & (~results[condition]["metrics"]["TPA"].isnull())
            ]
            # classnames.append(classname)
            results[condition]["metrics"].loc[
                results_class.index, "CA_relevant"
            ] = "yes"
            results[condition]["class_abundance"][classname] = {}
            results[condition]["class_abundance"][classname]["CA"] = np.median(
                results_class["TPA"]
            )
            results[condition]["class_abundance"][classname]["std"] = np.std(
                results_class["TPA"]
            )
            results[condition]["class_abundance"][classname]["count"] = len(
                results_class
            )

        # add CC:

        stacked_CC = np.stack(learning_xyz[condition]["z_full"].values())
        stacked_CC = np.array(
            [array / array.sum(axis=1, keepdims=True) for array in stacked_CC]
        )

        median_CC = np.median(stacked_CC, axis=0)
        Q1_CC = np.percentile(stacked_CC, 25, axis=0)
        Q3_CC = np.percentile(stacked_CC, 75, axis=0)

        normalization_factor = median_CC.sum(axis=1, keepdims=True)

        normalized_median_CC = median_CC / normalization_factor
        normalized_Q1_CC = Q1_CC / normalization_factor
        normalized_Q3_CC = Q3_CC / normalization_factor

        IQR_CC = normalized_Q3_CC - normalized_Q1_CC

        median_CC_df = pd.DataFrame(
            normalized_median_CC,
            columns=list(learning_xyz[condition]["z_full_df"].values())[
                0
            ].columns,
            index=list(learning_xyz[condition]["z_full_df"].values())[0].index,
        ).add_prefix("CC_")
        IQR_CC_df = pd.DataFrame(
            IQR_CC,
            columns=list(learning_xyz[condition]["z_full_df"].values())[
                0
            ].columns,
            index=list(learning_xyz[condition]["z_full_df"].values())[0].index,
        ).add_prefix("dCC_")

        results[condition]["metrics"] = pd.merge(
            results[condition]["metrics"],
            median_CC_df,
            left_index=True,
            right_index=True,
            how="left",
        )
        results[condition]["metrics"] = pd.merge(
            results[condition]["metrics"],
            IQR_CC_df,
            left_index=True,
            right_index=True,
            how="left",
        )

        for class_act in learning_xyz[condition]["classes"]:
            nonmarker_z = results[condition]["metrics"].loc[
                (results[condition]["metrics"]["marker"] != class_act)
                & (results[condition]["metrics"]["marker"].isna() == False)
            ][["CC_" + class_act]]
            thresh = np.percentile(
                nonmarker_z["CC_" + class_act].tolist(),
                NN_params["reliability"],
            )
            results[condition]["metrics"]["fCC_" + class_act] = results[
                condition
            ]["metrics"]["CC_" + class_act]
            results[condition]["metrics"].loc[
                results[condition]["metrics"]["fCC_" + class_act] < thresh,
                "fCC_" + class_act,
            ] = 0.0

        fCC_columns = [
            col
            for col in results[condition]["metrics"].columns
            if col.startswith("fCC_")
        ]
        row_sums = results[condition]["metrics"][fCC_columns].sum(axis=1)
        results[condition]["metrics"][fCC_columns] = results[condition][
            "metrics"
        ][fCC_columns].div(row_sums, axis=0)

        # add nCC:
        # classnames = list(learning_xyz[condition]['x_full_df'].columns)
        classnames = list(learning_xyz[condition]["z_full_df"].values())[
            0
        ].columns
        num_rows = stacked_CC.shape[1]

        # Create a dictionary to hold the data
        data = {"CClist_" + classname: [] for classname in classnames}

        # Populate the dictionary with lists for each row and column
        for col_idx, classname in enumerate(classnames):
            for row_idx in range(num_rows):
                row_values = stacked_CC[:, row_idx, col_idx]

                # Check if the entire row is NaN and assign NaN or the list of values
                if np.isnan(row_values).all():
                    data["CClist" + classname].append(np.nan)
                else:
                    data["CClist_" + classname].append(row_values.tolist())

        # Create a DataFrame from the dictionary
        list_df = pd.DataFrame(
            data,
            index=list(learning_xyz[condition]["z_full_df"].values())[0].index,
        )

        results[condition]["metrics"] = pd.merge(
            results[condition]["metrics"],
            list_df,
            left_index=True,
            right_index=True,
            how="left",
        )
        results[condition]["classnames"] = classnames

        for classname in classnames:
            results[condition]["metrics"]["nCC_" + classname] = (
                results[condition]["metrics"]["fCC_" + classname]
                * results[condition]["class_abundance"][classname]["CA"]
            )

        nCC_cols = [
            col
            for col in results[condition]["metrics"]
            if col.startswith("nCC_")
        ]
        nCC_sums = results[condition]["metrics"][nCC_cols].sum(axis=1)
        nCC_sums[nCC_sums == 0] = 1
        results[condition]["metrics"][nCC_cols] = results[condition][
            "metrics"
        ][nCC_cols].div(nCC_sums, axis=0)

        for classname in classnames:
            results[condition]["metrics"]["CPA_" + classname] = (
                results[condition]["metrics"]["CC_" + classname]
                * results[condition]["metrics"]["TPA"]
            )

            results[condition]["metrics"]["nCPA_" + classname] = (
                results[condition]["metrics"]["nCC_" + classname]
                * results[condition]["metrics"]["TPA"]
            )

        # refix = 'your_prefix_'
        # prefix_cols = [col for col in df.columns if col.startswith(prefix)]

        # # Calculate row-wise sums for the prefixed columns
        # row_sums = df[prefix_cols].sum(axis=1)

        # # Avoid division by zero
        # row_sums[row_sums == 0] = 1

        # # Divide each of the prefixed columns by the row-wise sums
        # df[prefix_cols] = df[prefix_cols].div(row_sums, axis=0)

        # add CC:
        # Assume your dictionary of NumPy arrays is named 'array_dict'
    # Example: array_dict = {'array1': array1, 'array2': array2, ...}

    # Stack the arrays along a new dimension
    # stacked_arrays = np.stack(array_dict.values())

    # # Calculate median and standard deviation along the new dimension (axis=0)
    # median_array = np.median(stacked_arrays, axis=0)
    # std_dev_array = np.std(stacked_arrays, axis=0)

    # # If you need the result back in DataFrame format:
    # median_df = pd.DataFrame(median_array, columns=list(df_dict.values())[0].columns, index=list(df_dict.values())[0].index)
    # std_dev_df = pd.DataFrame(std_dev_array, columns=list(df_dict.values())[0].columns, index=list(df_dict.values())[0].index)

    # print(results)

    # print('calculate TPA...')
    # TPA_list = []
    # for replicate in tp_data[condition]:
    #     TPA_list
    return results


def compare_lists(list1, list2):
    # Function to compare two lists and handle NaNs
    return sum(
        abs(a - b)
        for a, b in zip(list1, list2)
        if not pd.isna(a) and not pd.isna(b)
    )


def comp_exec3(mode, results, learning_xyz):
    conditions = [x for x in results]
    for condition in conditions:
        results[condition]["metrics"] = results[condition]["metrics"][
            ~results[condition]["metrics"].index.duplicated(keep="first")
        ]

    # class_lists = []
    # for condition in conditions:
    #     class_lists.append(results[condition]['classnames'])
    # classnames = list(set(class_lists[0]).intersection(*class_lists[1:]))

    combinations = []
    for con_1 in conditions:
        for con_2 in conditions:
            if con_1 != con_2:
                combinations.append((con_1, con_2))

    comparison = {}
    for comb in combinations:
        print(comb)
        classnames = list(
            set(results[comb[0]]["classnames"])
            & set(results[comb[1]]["classnames"])
        )
        comparison[comb] = {}

        ## prepare data:
        comparison[comb]["intersection_data"] = pd.merge(
            results[comb[0]]["metrics"],
            results[comb[1]]["metrics"],
            left_index=True,
            right_index=True,
            how="inner",
        )
        comparison[comb]["metrics"] = pd.DataFrame(
            index=comparison[comb]["intersection_data"].index
        )

        metrics_own = results[comb[0]]["metrics"]
        metrics_other = results[comb[1]]["metrics"]

        print("performing t-tests..")

        ## create RL, nRL, RLS, and nRLS:
        for classname in classnames:
            comparison[comb]["metrics"]["RL_" + classname] = (
                results[comb[1]]["metrics"]["CC_" + classname]
                - results[comb[0]]["metrics"]["CC_" + classname]
            )
            if mode == "deep":
                comparison[comb]["metrics"]["nRL_" + classname] = (
                    results[comb[1]]["metrics"]["nCC_" + classname]
                    - results[comb[0]]["metrics"]["nCC_" + classname]
                )
            elif mode == "rough":
                pass
        rl_cols = [
            col
            for col in comparison[comb]["metrics"].columns
            if col.startswith("RL_")
        ]
        comparison[comb]["metrics"]["RLS"] = (
            comparison[comb]["metrics"][rl_cols].abs().sum(axis=1)
        )
        if mode == "deep":
            nrl_cols = [
                col
                for col in comparison[comb]["metrics"].columns
                if col.startswith("nRL_")
            ]
            comparison[comb]["metrics"]["nRLS"] = (
                comparison[comb]["metrics"][nrl_cols].abs().sum(axis=1)
            )
        elif mode == "rough":
            pass

        for classname in classnames:
            comparison[comb]["metrics"]["fRL_" + classname] = (
                results[comb[1]]["metrics"]["fCC_" + classname]
                - results[comb[0]]["metrics"]["fCC_" + classname]
            )
        frl_cols = [
            col
            for col in comparison[comb]["metrics"].columns
            if col.startswith("fRL_")
        ]
        comparison[comb]["metrics"]["fRLS"] = (
            comparison[comb]["metrics"][frl_cols].abs().sum(axis=1)
        )

        ## t-test for each class (with CClist)
        # and Cohen's D for each class
        # t_df, common_indices, ncc_columns = perform_t_tests_per_cell(metrics_own, metrics_other, 'CClist_')
        # for classname in classnames:
        #     t_df.rename(columns = {'CClist_'+classname+'_P' : 'P_'+classname,
        #                            'CClist_'+classname+'_T' : 'T_'+classname,
        #                            'CClist_'+classname+'_D' : 'D_'+classname}, inplace = True)

        #     # calculate DS:
        # d_columns = [col for col in t_df.columns if col.startswith('D_')]
        # t_df['DS'] = t_df[d_columns].abs().sum(axis=1)

        #     # add statistics to metrics:
        # comparison[comb]['metrics'] = pd.merge(comparison[comb]['metrics'], t_df, left_index = True, right_index = True, how = 'left')

        # u_df, common_indices, ncc_columns = perform_mann_whitney_tests_per_cell(metrics_own, metrics_other, 'CClist_')
        # t_df, common_indices, ncc_columns = perform_t_tests_per_cell(metrics_own, metrics_other, 'CClist_')

        test_df, common_indices, ncc_columns = (
            perform_mann_whitney_t_tests_per_cell(
                metrics_own, metrics_other, "CClist_"
            )
        )
        # for classname in classnames:
        #     u_df.rename(columns = {'CClist_'+classname+'_P' : 'P(u)_'+classname,
        #                            'CClist_'+classname+'_U' : 'U(u)_'+classname,
        #                            'CClist_'+classname+'_D' : 'D(u)_'+classname}, inplace = True)
        for classname in classnames:
            test_df.rename(
                columns={
                    "CClist_" + classname + "_U": "U_" + classname,
                    "CClist_" + classname + "_T": "T_" + classname,
                    "CClist_" + classname + "_D": "D_" + classname,
                    "CClist_" + classname + "_P(U)": "P(U)_" + classname,
                    "CClist_" + classname + "_P(T)": "P(T)_" + classname,
                },
                inplace=True,
            )
            # calculate DS:
        d_columns = [col for col in test_df.columns if col.startswith("D_")]
        test_df["DS"] = test_df[d_columns].abs().sum(axis=1)

        # add statistics to metrics:
        comparison[comb]["metrics"] = pd.merge(
            comparison[comb]["metrics"],
            test_df,
            left_index=True,
            right_index=True,
            how="left",
        )

        print("calculate RLS lists...")
        RLS_results = {}
        RLS_null = {}
        for ID in common_indices:
            cclists_own = [
                metrics_own.loc[ID, "CClist_" + classname]
                for classname in classnames
            ]
            cclists_other = [
                metrics_other.loc[ID, "CClist_" + classname]
                for classname in classnames
            ]

            cclists_own_transposed = [
                list(values) for values in zip(*cclists_own)
            ]
            cclists_other_transposed = [
                list(values) for values in zip(*cclists_other)
            ]

            RLS_results[ID] = []
            RLS_null[ID] = []

            for i in range(len(cclists_own_transposed)):
                for j in range(i + 1, len(cclists_own_transposed)):
                    null_result = compare_lists(
                        cclists_own_transposed[i], cclists_own_transposed[j]
                    )
                    RLS_null[ID].append(null_result)
            for i in range(len(cclists_other_transposed)):
                for j in range(i + 1, len(cclists_other_transposed)):
                    null_result = compare_lists(
                        cclists_other_transposed[i],
                        cclists_other_transposed[j],
                    )
                    RLS_null[ID].append(null_result)

            for own_list in cclists_own_transposed:
                for other_list in cclists_other_transposed:
                    comparison_result = compare_lists(own_list, other_list)
                    RLS_results[ID].append(comparison_result)
        comparison[comb]["RLS_results"] = pd.Series(RLS_results)
        comparison[comb]["RLS_null"] = pd.Series(RLS_null)

        comparison[comb]["metrics"]["P(t)_RLS"] = np.nan
        comparison[comb]["metrics"]["P(u)_RLS"] = np.nan
        for index in comparison[comb]["metrics"].index:
            if index in common_indices:
                # Perform the t-test
                stat, p_value = ttest_ind(
                    comparison[comb]["RLS_results"].loc[index],
                    comparison[comb]["RLS_null"].loc[index],
                    nan_policy="omit",
                )
                comparison[comb]["metrics"].loc[index, "P(t)_RLS"] = p_value
                if (
                    is_all_nan(comparison[comb]["RLS_results"].loc[index])
                    or is_all_nan(comparison[comb]["RLS_null"].loc[index])
                    or len(set(comparison[comb]["RLS_results"].loc[index]))
                    == 1
                    or len(set(comparison[comb]["RLS_null"].loc[index])) == 1
                ):
                    comparison[comb]["metrics"].loc[index, "P(u)_RLS"] = pd.NA
                else:
                    stat_u, p_value_u = stats.mannwhitneyu(
                        comparison[comb]["RLS_results"].loc[index],
                        comparison[comb]["RLS_null"].loc[index],
                        alternative="two-sided",
                    )
                    comparison[comb]["metrics"].loc[index, "P(u)_RLS"] = (
                        p_value_u
                    )
            else:
                comparison[comb]["metrics"].loc[index, "P(t)_RLS"] = pd.NA
                comparison[comb]["metrics"].loc[index, "P(u)_RLS"] = pd.NA

        print("calculate nRLS lists...")
        if mode == "deep":
            nRLS_results = {}
            nRLS_null = {}
            for ID in common_indices:
                ncclists_own = [
                    metrics_own.loc[ID, "nCClist_" + classname]
                    for classname in classnames
                ]
                ncclists_other = [
                    metrics_other.loc[ID, "nCClist_" + classname]
                    for classname in classnames
                ]

                ncclists_own_transposed = [
                    list(values) for values in zip(*ncclists_own)
                ]
                ncclists_other_transposed = [
                    list(values) for values in zip(*ncclists_other)
                ]

                nRLS_results[ID] = []
                nRLS_null[ID] = []

                for i in range(len(ncclists_own_transposed)):
                    for j in range(i + 1, len(ncclists_own_transposed)):
                        null_result = compare_lists(
                            ncclists_own_transposed[i],
                            ncclists_own_transposed[j],
                        )
                        nRLS_null[ID].append(null_result)
                for i in range(len(ncclists_other_transposed)):
                    for j in range(i + 1, len(ncclists_other_transposed)):
                        null_result = compare_lists(
                            ncclists_other_transposed[i],
                            ncclists_other_transposed[j],
                        )
                        nRLS_null[ID].append(null_result)

                for own_list in ncclists_own_transposed:
                    for other_list in ncclists_other_transposed:
                        comparison_result = compare_lists(own_list, other_list)
                        nRLS_results[ID].append(comparison_result)
            comparison[comb]["nRLS_results"] = pd.Series(nRLS_results)
            comparison[comb]["nRLS_null"] = pd.Series(nRLS_null)

            comparison[comb]["metrics"]["P(t)_nRLS"] = np.nan
            comparison[comb]["metrics"]["P(u)_nRLS"] = np.nan
            for index in comparison[comb]["metrics"].index:
                if index in common_indices:
                    # Perform the t-test
                    stat, p_value = ttest_ind(
                        comparison[comb]["nRLS_results"].loc[index],
                        comparison[comb]["nRLS_null"].loc[index],
                        nan_policy="omit",
                    )
                    comparison[comb]["metrics"].loc[index, "P(t)_nRLS"] = (
                        p_value
                    )
                    if (
                        is_all_nan(comparison[comb]["nRLS_results"].loc[index])
                        or is_all_nan(comparison[comb]["nRLS_null"].loc[index])
                        or len(
                            set(comparison[comb]["nRLS_results"].loc[index])
                        )
                        == 1
                        or len(set(comparison[comb]["nRLS_null"].loc[index]))
                        == 1
                    ):
                        comparison[comb]["metrics"].loc[index, "P(u)_nRLS"] = (
                            pd.NA
                        )
                    else:
                        stat_u, p_value_u = stats.mannwhitneyu(
                            comparison[comb]["nRLS_results"].loc[index],
                            comparison[comb]["nRLS_null"].loc[index],
                            alternative="two-sided",
                        )
                        comparison[comb]["metrics"].loc[index, "P(u)_nRLS"] = (
                            p_value_u
                        )
                else:
                    comparison[comb]["metrics"].loc[index, "P(t)_nRLS"] = pd.NA
                    comparison[comb]["metrics"].loc[index, "P(u)_nRLS"] = pd.NA
        elif mode == "rough":
            pass

        print("calculate CPAs...")
        if mode == "deep":
            for classname in classnames:
                metrics_own["CPA_log_" + classname] = np.log2(
                    metrics_own["CPA_" + classname]
                )
                metrics_own = impute_data(
                    metrics_own, "CPA_log_" + classname, "CPA_imp_" + classname
                )

                metrics_other["CPA_log_" + classname] = np.log2(
                    metrics_other["CPA_" + classname]
                )
                metrics_other = impute_data(
                    metrics_other,
                    "CPA_log_" + classname,
                    "CPA_imp_" + classname,
                )

                comparison[comb]["metrics"]["CFC_" + classname] = (
                    metrics_other["CPA_imp_" + classname]
                    - metrics_own["CPA_imp_" + classname]
                )

                metrics_own["nCPA_log_" + classname] = np.log2(
                    metrics_own["nCPA_" + classname]
                )
                metrics_own = impute_data(
                    metrics_own,
                    "nCPA_log_" + classname,
                    "nCPA_imp_" + classname,
                )

                metrics_other["nCPA_log_" + classname] = np.log2(
                    metrics_other["nCPA_" + classname]
                )
                metrics_other = impute_data(
                    metrics_other,
                    "nCPA_log_" + classname,
                    "nCPA_imp_" + classname,
                )

                comparison[comb]["metrics"]["nCFC_" + classname] = (
                    metrics_other["nCPA_imp_" + classname]
                    - metrics_own["nCPA_imp_" + classname]
                )
        elif mode == "rough":
            pass

    return comparison


def comp_exec2(results, learning_xyz):
    conditions = [x for x in results]
    for condition in conditions:
        results[condition]["metrics"] = results[condition]["metrics"][
            ~results[condition]["metrics"].index.duplicated(keep="first")
        ]

    # class_lists = []
    # for condition in conditions:
    #     class_lists.append(results[condition]['classnames'])
    # classnames = list(set(class_lists[0]).intersection(*class_lists[1:]))

    combinations = []
    for con_1 in conditions:
        for con_2 in conditions:
            if con_1 != con_2:
                combinations.append((con_1, con_2))

    comparison = {}
    for comb in combinations:
        print(comb)
        # classnames = list(set(results[comb[0]]['classnames'] & results[comb[1]]['classnames']))
        classnames = list(
            set(results[comb[0]]["classnames"])
            & set(results[comb[1]]["classnames"])
        )
        comparison[comb] = {}

        # subcons_own = [x for x in learning_xyz if x.startswith(comb[0]+'_')]
        # subcons_other = [x for x in learning_xyz if x.startswith(comb[1]+'_')]
        # x_full_df_own = pd.DataFrame()
        # for subcon in subcons_own:
        #     x_full_df_own = pd.merge(x_full_df_own, learning_xyz[subcon]['x_full_df'], left_index = True, right_index = True, how = 'inner')
        # x_full_df_other = pd.DataFrame()
        # for subcon in subcons_other:
        #     x_

        comparison[comb]["intersection_data"] = pd.merge(
            results[comb[0]]["metrics"],
            results[comb[1]]["metrics"],
            left_index=True,
            right_index=True,
            how="inner",
        )
        comparison[comb]["metrics"] = pd.DataFrame(
            index=comparison[comb]["intersection_data"].index
        )

        # comparison[comb]['intersection_data'] = pd.merge(learning_xyz[comb[0]]['x_full_df'], learning_xyz[comb[1]]['x_full_df'], left_index = True, right_index = True, how = 'inner')
        # comparison[comb]['metrics'] = pd.DataFrame(index = comparison[comb]['intersection_data'].index)
        # comparison[comb]['nCC'] = pd.DataFrame(index = comparison[comb]['intersection_data'].index)

        metrics_own = results[comb[0]]["metrics"]
        metrics_other = results[comb[1]]["metrics"]

        print("performing t-tests..")
        ##### for nCC:
        # paired_t_df, common_indices, ncc_columns = perform_paired_t_tests_per_cell(metrics_own, metrics_other, 'CClist_')
        # comparison[comb]['nCC'] = pd.merge(comparison[comb]['nCC'], paired_t_df, left_index = True, right_index = True, how = 'left')

        for classname in classnames:
            comparison[comb]["metrics"]["RL_" + classname] = (
                results[comb[1]]["metrics"]["CC_" + classname]
                - results[comb[0]]["metrics"]["CC_" + classname]
            )
            comparison[comb]["metrics"]["nRL_" + classname] = (
                results[comb[1]]["metrics"]["nCC_" + classname]
                - results[comb[0]]["metrics"]["nCC_" + classname]
            )
        rl_cols = [
            col
            for col in comparison[comb]["metrics"].columns
            if col.startswith("RL_")
        ]
        nrl_cols = [
            col
            for col in comparison[comb]["metrics"].columns
            if col.startswith("nRL_")
        ]
        comparison[comb]["metrics"]["RLS"] = (
            comparison[comb]["metrics"][rl_cols].abs().sum(axis=1)
        )
        comparison[comb]["metrics"]["nRLS"] = (
            comparison[comb]["metrics"][nrl_cols].abs().sum(axis=1)
        )

        ##### for CC:
        # paired_t_df, common_indices, ncc_columns = perform_paired_t_tests_per_cell(metrics_own, metrics_other, 'CClist_')
        t_df, common_indices, ncc_columns = perform_t_tests_per_cell(
            metrics_own, metrics_other, "CClist_"
        )
        # t_stat, p_value = ttest_ind(df_1['col_n'].dropna(), df_2['col_n'].dropna(), equal_var=False)

        for classname in classnames:
            t_df.rename(
                columns={
                    "CClist_" + classname + "_P": "P_" + classname,
                    "CClist_" + classname + "_T": "T_" + classname,
                    "CClist_" + classname + "_D": "D_" + classname,
                },
                inplace=True,
            )

        d_columns = [col for col in t_df.columns if col.startswith("D_")]
        # paired_t_df['RLS'] = paired_t_df[d_columns].sum(axis=1)
        t_df["DS"] = t_df[d_columns].abs().sum(axis=1)
        comparison[comb]["metrics"] = pd.merge(
            comparison[comb]["metrics"],
            t_df,
            left_index=True,
            right_index=True,
            how="left",
        )

        for classname in classnames:
            col_old = "P_" + classname
            col_new = "FDR_" + classname

            comparison[comb]["metrics"] = benjamini_hochberg_correction(
                comparison[comb]["metrics"], col_old, col_new
            )

            # colname = 'CC_'+classname+'_P'
            # comparison[comb]['metrics']['CC_'+classname+'_FDR'] = benjamini_hochberg_correction(comparison[comb]['metrics'], colname)

        print("performing hotelling tests...")

        ##### for CC:

        # t_statistics, p_statistics = process_and_test(metrics_own, metrics_other, classnames)
        t2_results = process_and_test(metrics_own, metrics_other, classnames)

        for idx, data in t2_results.items():
            comparison[comb]["metrics"].at[idx, "t2_T"] = data["t_stat"]
            comparison[comb]["metrics"].at[idx, "t2_P"] = data["p_value"]
        # comparison[comb]['metrics']['t2_T'] = np.nan
        # comparison[comb]['metrics']['t2_P'] = np.nan
        # for idx in common_indices:
        #     row_df1_data = []
        #     row_df2_data = []
        #     for col in ncc_columns:
        #         if isinstance(metrics_own.at[idx, col], list) and isinstance(metrics_other.at[idx, col], list):
        #             if metrics_own.at[idx, col] and metrics_other.at[idx, col] and not all(np.isnan(metrics_own.at[idx, col])) and not all(np.isnan(metrics_other.at[idx, col])):
        #                 imputed_data_df1 = impute_nan_with_mean(metrics_own.at[idx, col])
        #                 imputed_data_df2 = impute_nan_with_mean(metrics_other.at[idx, col])
        #                 row_df1_data.append(imputed_data_df1)
        #                 row_df2_data.append(imputed_data_df2)

        #     if row_df1_data and row_df2_data:
        #         row_df1 = np.array(row_df1_data).T
        #         row_df2 = np.array(row_df2_data).T
        #         if row_df1.shape == row_df2.shape and row_df1.size > 0 and row_df2.size > 0:
        #             #t2_stat, t2_p_value = hotelling_t2_test(row_df1, row_df2)
        #             t2_stat, t2_p_value, t2_diff = hotelling_t2_test_rowwise(row_df1, row_df2)

        #             comparison[comb]['metrics'].at[idx, 't2_T'] = t2_stat
        #             comparison[comb]['metrics'].at[idx, 't2_P'] = t2_p_value
        #         else:
        #             comparison[comb]['metrics'].at[idx, 't2_T'] = np.nan
        #             comparison[comb]['metrics'].at[idx, 't2_P'] = np.nan
        #     else:
        #         comparison[comb]['metrics'].at[idx, 't2_T'] = np.nan
        #         comparison[comb]['metrics'].at[idx, 't2_P'] = np.nan

        comparison[comb]["metrics"] = benjamini_hochberg_correction(
            comparison[comb]["metrics"], "t2_P", "t2_FDR"
        )

        for classname in classnames:
            # CPA_col = 'CPA_log_'+classname
            metrics_own["CPA_log_" + classname] = np.log2(
                metrics_own["CPA_" + classname]
            )
            metrics_own = impute_data(
                metrics_own, "CPA_log_" + classname, "CPA_imp_" + classname
            )

            metrics_own["nCPA_log_" + classname] = np.log2(
                metrics_own["nCPA_" + classname]
            )
            metrics_own = impute_data(
                metrics_own, "nCPA_log_" + classname, "nCPA_imp_" + classname
            )

            metrics_other["CPA_log_" + classname] = np.log2(
                metrics_other["CPA_" + classname]
            )
            metrics_other = impute_data(
                metrics_other, "CPA_log_" + classname, "CPA_imp_" + classname
            )

            metrics_other["nCPA_log_" + classname] = np.log2(
                metrics_other["nCPA_" + classname]
            )
            metrics_other = impute_data(
                metrics_other, "nCPA_log_" + classname, "nCPA_imp_" + classname
            )

            # nCPA_col = 'nCPA_log_'+classname
            # metrics_own['nCPA_log_'+classname] = np.log2(metrics_own['nCPA_'+classname])
            # metrics_other = impute_data(metrics_other, 'CPA_log_'+classname, 'CPA_imp_'+classname)

            # impute_data (results[condition]['metrics'], 'CPA_'+classname)

            comparison[comb]["metrics"]["CFC_" + classname] = (
                metrics_other["CPA_imp_" + classname]
                - metrics_own["CPA_imp_" + classname]
            )
            comparison[comb]["metrics"]["nCFC_" + classname] = (
                metrics_other["nCPA_imp_" + classname]
                - metrics_own["nCPA_imp_" + classname]
            )

            # comparison[comb]['metrics']['CFC'+classname] = np.log2(metrics_other['CPA_'+classname]) - np.log2(metrics_own['CPA_'+classname])
            # comparison[comb]['metrics']['nCFC'+classname] = np.log2(metrics_other['nCPA_'+classname]) - np.log2(metrics_own['nCPA_'+classname])

    return comparison


# def hotelling_t2_test(arr1, arr2, regularization_value=1e-6):
#     if arr1.size == 0 or arr2.size == 0:
#         return np.nan, np.nan

#     # Calculate means
#     mean1 = np.nanmean(arr1, axis=0)
#     mean2 = np.nanmean(arr2, axis=0)
#     mean_diff = mean1 - mean2

#     # Calculate pooled covariance matrix with regularization
#     n1, n2 = arr1.shape[0], arr2.shape[0]
#     cov1 = np.nan_to_num(np.nanvar(arr1, axis=0, ddof=1))  # Convert NaNs to zero
#     cov2 = np.nan_to_num(np.nanvar(arr2, axis=0, ddof=1))  # Convert NaNs to zero
#     pooled_cov = ((n1 - 1) * cov1 + (n2 - 1) * cov2) / (n1 + n2 - 2)
#     pooled_cov += np.eye(pooled_cov.shape[0]) * regularization_value

#     if np.linalg.cond(pooled_cov) > 1 / np.finfo(pooled_cov.dtype).eps:
#         return np.nan, np.nan  # Skip if the matrix is ill-conditioned

#     t_squared = n1 * n2 / (n1 + n2) * np.dot(mean_diff.T, np.dot(np.linalg.pinv(pooled_cov), mean_diff))

#     p = arr1.shape[1]  # number of variables
#     df1 = p
#     df2 = n1 + n2 - p - 1

#     if df1 <= 0 or df2 <= 0:
#         return np.nan, np.nan

#     f_stat = (t_squared / p) * (df2 / df1)
#     p_value = 1 - f.cdf(f_stat, df1, df2)

#     return t_squared, p_value


# def hotelling_t2_test(arr1, arr2, regularization_value=1e-6):
#     # Ensure arr1 and arr2 are two-dimensional
#     arr1 = np.atleast_2d(arr1)
#     arr2 = np.atleast_2d(arr2)

#     if arr1.size == 0 or arr2.size == 0 or arr1.shape[1] != arr2.shape[1]:
#         # If either array is empty or they don't have the same number of columns
#         return np.nan, np.nan

#     n1, n2 = arr1.shape[0], arr2.shape[0]

#     # Calculate means, considering NaNs
#     mean1 = np.nanmean(arr1, axis=0)
#     mean2 = np.nanmean(arr2, axis=0)
#     mean_diff = mean1 - mean2

#     # Covariance matrices
#     cov1 = np.nanvar(arr1, axis=0, ddof=1)
#     cov2 = np.nanvar(arr2, axis=0, ddof=1)

#     # Adjust shapes for covariance matrices if they are not 2D
#     cov1 = np.atleast_2d(cov1)
#     cov2 = np.atleast_2d(cov2)

#     # Calculate pooled covariance
#     pooled_cov = ((n1 - 1) * cov1 + (n2 - 1) * cov2) / (n1 + n2 - 2)

#     # Regularize pooled covariance
#     if pooled_cov.shape == (1, 1):
#         pooled_cov += regularization_value
#     elif pooled_cov.shape[0] > 1:
#         pooled_cov += np.eye(pooled_cov.shape[0]) * regularization_value
#     else:
#         # If pooled_cov is not a square matrix, return NaN
#         return np.nan, np.nan

#     # Calculate inverse of pooled covariance matrix
#     try:
#         inverse_pooled_cov = np.linalg.pinv(pooled_cov)
#     except np.linalg.LinAlgError:
#         return np.nan, np.nan

#     # Hotelling's T-squared statistic
#     t_squared = n1 * n2 / (n1 + n2) * mean_diff.T @ inverse_pooled_cov @ mean_diff

#     # Degrees of freedom
#     p = arr1.shape[1]  # Number of variables
#     df1 = p
#     df2 = n1 + n2 - p - 1

#     if df1 <= 0 or df2 <= 0:
#         return np.nan, np.nan

#     # Calculate F-statistic and p-value
#     f_stat = (t_squared / p) * (df2 / df1)
#     p_value = 1 - f.cdf(f_stat, df1, df2)


#     return t_squared, p_value
def hotelling_t2_test(data1, data2, regularization_value=1e-6):
    mean1 = np.mean(data1, axis=0)
    mean2 = np.mean(data2, axis=0)

    n1, n2 = data1.shape[0], data2.shape[0]
    cov1 = np.cov(data1, rowvar=False)
    cov2 = np.cov(data2, rowvar=False)

    # Ensure pooled_cov is 2D
    pooled_cov = np.atleast_2d(
        ((n1 - 1) * cov1 + (n2 - 1) * cov2) / (n1 + n2 - 2)
    )
    pooled_cov += np.eye(pooled_cov.shape[0]) * regularization_value

    mean_diff = mean1 - mean2
    t_squared = (
        n1
        * n2
        / (n1 + n2)
        * np.dot(mean_diff.T, np.dot(np.linalg.inv(pooled_cov), mean_diff))
    )

    p = data1.shape[1]  # number of variables
    df1 = p
    df2 = n1 + n2 - p - 1

    f_stat = (df2 / (n1 + n2 - p - 1)) * t_squared
    p_value = 1 - f.cdf(f_stat, df1, df2)

    return t_squared, p_value


def process_and_test(data_own, data_other, classnames):
    results = {}

    filtered_cols = [
        col
        for col in data_own.columns
        if col.startswith("CClist_") and col.split("_")[1] in classnames
    ]

    for idx in data_own.index:
        # Check if index exists in both dataframes
        if idx not in data_other.index:
            results[idx] = {"t_stat": np.nan, "p_value": np.nan}
            continue

        valid_cols = [
            col
            for col in filtered_cols
            if col in data_own.columns and col in data_other.columns
        ]
        row_df1_data, row_df2_data = [], []

        for col in valid_cols:
            val1, val2 = data_own.at[idx, col], data_other.at[idx, col]
            if (
                isinstance(val1, list)
                and isinstance(val2, list)
                and not all(np.isnan(val1))
                and not all(np.isnan(val2))
            ):
                row_df1_data.extend(val1)
                row_df2_data.extend(val2)

        if len(row_df1_data) > 0 and len(row_df2_data) > 0:
            row_df1 = np.array(row_df1_data).reshape(-1, 1)
            row_df2 = np.array(row_df2_data).reshape(-1, 1)
            t_stat, p_val = hotelling_t2_test(row_df1, row_df2)
            results[idx] = {"t_stat": t_stat, "p_value": p_val}
        else:
            results[idx] = {"t_stat": np.nan, "p_value": np.nan}

    return results


# def hotelling_t2_test_rowwise(data1, data2, regularization_value=1e-6):
#     t_squared_results = []
#     p_value_results = []
#     mean_diff_results = []

#     for row1, row2 in zip(data1, data2):
#         # Flatten the lists, filter out NaNs, and convert to numpy arrays
#         arr1 = np.array([x for sublist in row1 for x in sublist if not np.isnan(x)])
#         arr2 = np.array([x for sublist in row2 for x in sublist if not np.isnan(x)])

#         # Skip the row if either array is empty after removing NaNs
#         if arr1.size == 0 or arr2.size == 0:
#             t_squared_results.append(np.nan)
#             p_value_results.append(np.nan)
#             mean_diff_results.append(np.nan)
#             continue

#         # Calculate means
#         mean1 = np.mean(arr1, axis=0)
#         mean2 = np.mean(arr2, axis=0)
#         mean_diff = mean1 - mean2
#         mean_diff_results.append(mean_diff)

#         # Calculate the pooled covariance matrix
#         n1, n2 = arr1.shape[0], arr2.shape[0]
#         cov1 = np.cov(arr1, rowvar=False)
#         cov2 = np.cov(arr2, rowvar=False)
#         pooled_cov = ((n1 - 1) * cov1 + (n2 - 1) * cov2) / (n1 + n2 - 2)
#         pooled_cov += np.eye(pooled_cov.shape[0]) * regularization_value

#         # Calculate the Hotelling's T-squared statistic
#         t_squared = n1 * n2 / (n1 + n2) * np.dot(np.dot(mean_diff.T, np.linalg.inv(pooled_cov)), mean_diff)
#         t_squared_results.append(t_squared)

#         # Calculate the degrees of freedom for the distribution
#         p = arr1.shape[1]  # number of variables
#         df1 = p
#         df2 = n1 + n2 - p - 1

#         # Calculate the F-statistic and p-value
#         f_stat = (df2 / (n1 + n2 - p - 1)) * t_squared
#         p_value = 1 - f.cdf(f_stat, df1, df2)
#         p_value_results.append(p_value)

#     return t_squared_results, p_value_results, mean_diff_results


def impute_data(df, colname, newcol):
    s = 1.8
    w = 0.3

    mean = np.mean(
        df[colname].replace(-np.inf, np.nan)
    )  # Exclude 0s from mean calculation
    std = np.std(
        df[colname].replace(-np.inf, np.nan)
    )  # Exclude 0s from std calculation
    mean_imp = mean - s * std  # Use your specified 's' value
    sigma = std * w  # Use your specified 'w' value

    # Apply the imputation for 0 values
    df[newcol] = df[colname].apply(
        lambda x: np.random.normal(mean_imp, sigma, 1)[0]
        if x == -np.inf
        else x
    )

    return df


def benjamini_hochberg_correction(df, column_old, column_new):
    p_values = pd.to_numeric(df[column_old], errors="coerce")
    non_nan_p_values = p_values.dropna()
    sorted_indices = non_nan_p_values.argsort()
    sorted_p_values = non_nan_p_values.iloc[sorted_indices]

    total_tests = len(sorted_p_values)
    adjusted_p_values = np.full(df.shape[0], np.nan)  # Initialize with NaNs
    prev_adj_p_val = 0

    # Map pandas indices to range-based integer positions
    index_mapping = {index: i for i, index in enumerate(df.index)}

    for sorted_index, p_value in zip(sorted_indices.index, sorted_p_values):
        original_index = index_mapping[sorted_index]  # Get integer position
        adjusted_p = min(p_value * total_tests / (original_index + 1), 1)
        adjusted_p = max(adjusted_p, prev_adj_p_val)
        prev_adj_p_val = adjusted_p
        adjusted_p_values[original_index] = adjusted_p

    adjusted_p_values_series = pd.Series(adjusted_p_values, index=df.index)
    df[column_new] = adjusted_p_values_series

    return df


# def benjamini_hochberg_correction(df, column_old, column_new, fdr_threshold=0.05):
#     # Convert the column to numeric, coercing errors to NaN
#     p_values = pd.to_numeric(df[column_old], errors='coerce')

#     # Handling NaNs: Only consider non-NaN p-values for sorting and adjustment
#     non_nan_p_values = p_values.dropna()
#     sorted_indices = non_nan_p_values.argsort()
#     sorted_p_values = non_nan_p_values.iloc[sorted_indices]

#     total_tests = len(sorted_p_values)
#     adjusted_p_values = np.full(sorted_p_values.shape, np.nan)

#     for i, p_value in enumerate(sorted_p_values):
#         adjusted_p = p_value * total_tests / (i + 1)
#         adjusted_p_values[sorted_indices[i]] = adjusted_p

#     # Creating a new Series for adjusted p-values, aligning with the original DataFrame's index
#     adjusted_p_values_series = pd.Series(adjusted_p_values, index=sorted_indices.index)

#     # Assigning adjusted p-values to the DataFrame
#     df[column_new] = adjusted_p_values_series.reindex(df.index)

#     return df

# # Example usage:
# for comb in combinations:  # Assuming combinations is defined in your context
#     benjamini_hochberg_correction(comparison[comb]['comparison'], 't2_p-val')
#     # Now comparison[comb]['comparison'] contains adjusted p-values and a 'significant' flag


def impute_nan_with_mean(data):
    non_nan_data = [x for x in data if not np.isnan(x)]
    mean_value = np.mean(non_nan_data) if non_nan_data else np.nan
    return [x if not np.isnan(x) else mean_value for x in data]


# def hotelling_t2_test(data1, data2, regularization_value=1e-6):
#     # Calculate means
#     mean1 = np.mean(data1, axis=0)
#     mean2 = np.mean(data2, axis=0)

#     # Calculate the pooled covariance matrix
#     n1, n2 = data1.shape[0], data2.shape[0]
#     cov1 = np.cov(data1, rowvar=False)
#     cov2 = np.cov(data2, rowvar=False)
#     pooled_cov = ((n1 - 1) * cov1 + (n2 - 1) * cov2) / (n1 + n2 - 2)

#     # pooled_cov = ((n1 - 1) * cov1 + (n2 - 1) * cov2) / (n1 + n2 - 2)
#     pooled_cov += np.eye(pooled_cov.shape[0]) * regularization_value

#     # Calculate the Hotelling's T-squared statistic
#     mean_diff = mean1 - mean2
#     t_squared = n1 * n2 / (n1 + n2) * np.dot(np.dot(mean_diff.T, np.linalg.inv(pooled_cov)), mean_diff)

#     # Calculate the degrees of freedom for the distribution
#     p = data1.shape[1]  # number of variables
#     df1 = p
#     df2 = n1 + n2 - p - 1

#     # Calculate the F-statistic
#     f_stat = (df2 / (n1 + n2 - p - 1)) * t_squared
#     p_value = 1 - f.cdf(f_stat, df1, df2)

#     return t_squared, p_value


# def perform_paired_t_tests_per_cell(df1, df2, prefix):
#     # Find columns starting with 'nCC_' and present in both dataframes
#     cc_columns = [col for col in df1.columns if col.startswith(prefix) and col in df2.columns]

#     # Find common indices in both dataframes
#     common_indices = df1.index.intersection(df2.index)

#     # Prepare the DataFrame to store the results
#     result_columns = [f"{col}_T" for col in cc_columns] + [f"{col}_P" for col in cc_columns]
#     results_df = pd.DataFrame(index=common_indices, columns=result_columns)

#     # Perform paired t-tests for each 'nCC_' column and each common index
#     for col in cc_columns:
#         for idx in common_indices:
#             list_df1 = df1.loc[idx, col]
#             list_df2 = df2.loc[idx, col]

#             # Check if either list is nan or contains only nan values
#             if is_all_nan(list_df1) or is_all_nan(list_df2):
#                 results_df.loc[idx, f"{col}_T"] = np.nan
#                 results_df.loc[idx, f"{col}_P"] = np.nan
#             else:
#                 # Perform paired t-test on the lists in the cell
#                 t_stat, p_value = stats.ttest_rel(list_df1, list_df2)
#                 results_df.loc[idx, f"{col}_T"] = t_stat
#                 results_df.loc[idx, f"{col}_P"] = p_value

#     return results_df, common_indices, cc_columns


def is_all_nan(list_):
    return (
        all(np.isnan(x) for x in list_)
        if isinstance(list_, list)
        else np.isnan(list_)
    )


def perform_paired_t_tests_per_cell(df1, df2, prefix):
    cc_columns = [
        col
        for col in df1.columns
        if col.startswith(prefix) and col in df2.columns
    ]
    common_indices = df1.index.intersection(df2.index)

    # Adding columns for Cohen's d
    result_columns = (
        [f"{col}_T" for col in cc_columns]
        + [f"{col}_P" for col in cc_columns]
        + [f"{col}_D" for col in cc_columns]
    )
    results_df = pd.DataFrame(index=common_indices, columns=result_columns)

    for col in cc_columns:
        for idx in common_indices:
            list_df1 = df1.loc[idx, col]
            list_df2 = df2.loc[idx, col]

            if is_all_nan(list_df1) or is_all_nan(list_df2):
                results_df.loc[idx, f"{col}_T"] = np.nan
                results_df.loc[idx, f"{col}_P"] = np.nan
                results_df.loc[idx, f"{col}_D"] = np.nan
            else:
                # Calculate mean difference and standard deviation of differences
                diff = np.array(list_df1) - np.array(list_df2)
                mean_diff = np.mean(diff)
                std_diff = np.std(
                    diff, ddof=1
                )  # ddof=1 for sample standard deviation

                # Perform paired t-test
                t_stat, p_value = stats.ttest_rel(list_df1, list_df2)

                # Calculate Cohen's d
                cohen_d = mean_diff / std_diff if std_diff != 0 else np.nan

                # Storing results
                results_df.loc[idx, f"{col}_T"] = t_stat
                results_df.loc[idx, f"{col}_P"] = p_value
                results_df.loc[idx, f"{col}_D"] = cohen_d

    return results_df, common_indices, cc_columns


# def perform_mann_whitney_tests_per_cell(df1, df2, prefix):
#     cc_columns = [col for col in df1.columns if col.startswith(prefix) and col in df2.columns]
#     common_indices = df1.index.intersection(df2.index)

#     # Columns for Mann-Whitney U results and Cohen's d
#     result_columns = [f"{col}_U" for col in cc_columns] + [f"{col}_P" for col in cc_columns] + [f"{col}_D" for col in cc_columns]
#     results_df = pd.DataFrame(index=common_indices, columns=result_columns)

#     for col in cc_columns:
#         for idx in common_indices:
#             list_df1 = df1.loc[idx, col]
#             list_df2 = df2.loc[idx, col]

#             if is_all_nan(list_df1) or is_all_nan(list_df2) or len(set(list_df1)) == 1 or len(set(list_df2)) == 1:
#                 # Handling cases where one or both groups are all NaNs or have no variability
#                 results_df.loc[idx, f"{col}_U"] = np.nan
#                 results_df.loc[idx, f"{col}_P"] = np.nan
#                 results_df.loc[idx, f"{col}_D"] = np.nan
#             else:
#                 # Perform Mann-Whitney U test
#                 u_stat, p_value = stats.mannwhitneyu(list_df1, list_df2, alternative='two-sided')

#                 # Calculating Cohen's d
#                 diff = [value_1 - value_2 for value_1 in list_df1 for value_2 in list_df2]
#                 mean_diff = np.mean(diff)
#                 std_diff = np.std(diff, ddof=1)
#                 cohen_d = mean_diff / std_diff if std_diff != 0 else np.nan

#                 # Storing results
#                 results_df.loc[idx, f"{col}_U"] = u_stat
#                 results_df.loc[idx, f"{col}_P"] = p_value
#                 results_df.loc[idx, f"{col}_D"] = cohen_d

#     return results_df, common_indices, cc_columns


def perform_mann_whitney_t_tests_per_cell(df1, df2, prefix):
    cc_columns = [
        col
        for col in df1.columns
        if col.startswith(prefix) and col in df2.columns
    ]
    common_indices = df1.index.intersection(df2.index)

    # Columns for Mann-Whitney U results, t-test results, and Cohen's d
    result_columns = (
        [f"{col}_U" for col in cc_columns]
        + [f"{col}_P_U" for col in cc_columns]
        + [f"{col}_D" for col in cc_columns]
        + [f"{col}_T" for col in cc_columns]
        + [f"{col}_P_T" for col in cc_columns]
    )
    results_df = pd.DataFrame(index=common_indices, columns=result_columns)

    for col in cc_columns:
        for idx in common_indices:
            list_df1 = df1.loc[idx, col]
            list_df2 = df2.loc[idx, col]

            if (
                is_all_nan(list_df1)
                or is_all_nan(list_df2)
                or len(set(list_df1)) == 1
                or len(set(list_df2)) == 1
            ):
                # Handling cases where one or both groups are all NaNs or have no variability
                results_df.loc[idx, f"{col}_U"] = np.nan
                results_df.loc[idx, f"{col}_P_U"] = np.nan
                results_df.loc[idx, f"{col}_D"] = np.nan
                results_df.loc[idx, f"{col}_T"] = np.nan
                results_df.loc[idx, f"{col}_P_T"] = np.nan
            else:
                # Perform Mann-Whitney U test
                u_stat, p_value_u = stats.mannwhitneyu(
                    list_df1, list_df2, alternative="two-sided"
                )

                # Perform t-test
                t_stat, p_value_t = stats.ttest_ind(
                    list_df1, list_df2, equal_var=False, nan_policy="omit"
                )

                # Calculating Cohen's d
                diff = [
                    value_1 - value_2
                    for value_1 in list_df1
                    for value_2 in list_df2
                ]
                mean_diff = np.mean(diff)
                std_diff = np.std(diff, ddof=1)
                cohen_d = mean_diff / std_diff if std_diff != 0 else np.nan

                # Storing results
                results_df.loc[idx, f"{col}_U"] = u_stat
                results_df.loc[idx, f"{col}_P(U)"] = p_value_u
                results_df.loc[idx, f"{col}_D"] = cohen_d
                results_df.loc[idx, f"{col}_T"] = t_stat
                results_df.loc[idx, f"{col}_P(T)"] = p_value_t

    return results_df, common_indices, cc_columns


def perform_t_tests_per_cell(df1, df2, prefix):
    cc_columns = [
        col
        for col in df1.columns
        if col.startswith(prefix) and col in df2.columns
    ]
    common_indices = df1.index.intersection(df2.index)

    # Adding columns for Cohen's d
    result_columns = (
        [f"{col}_T" for col in cc_columns]
        + [f"{col}_P" for col in cc_columns]
        + [f"{col}_D" for col in cc_columns]
    )
    results_df = pd.DataFrame(index=common_indices, columns=result_columns)

    for col in cc_columns:
        for idx in common_indices:
            list_df1 = df1.loc[idx, col]
            list_df2 = df2.loc[idx, col]

            if is_all_nan(list_df1) or is_all_nan(list_df2):
                results_df.loc[idx, f"{col}_T"] = np.nan
                results_df.loc[idx, f"{col}_P"] = np.nan
                results_df.loc[idx, f"{col}_D"] = np.nan
            else:
                diff = []
                for value_1 in list_df1:
                    for value_2 in list_df2:
                        diff.append(value_1 - value_2)

                # Calculate mean difference and standard deviation of differences
                # diff = np.mean(list_df1) - np.mean(list_df2)
                mean_diff = np.mean(diff)
                std_diff = np.std(
                    diff, ddof=1
                )  # ddof=1 for sample standard deviation

                # Perform paired t-test
                t_stat, p_value = stats.ttest_ind(list_df1, list_df2)

                # Calculate Cohen's d
                cohen_d = mean_diff / std_diff if std_diff != 0 else np.nan

                # Storing results
                results_df.loc[idx, f"{col}_T"] = t_stat
                results_df.loc[idx, f"{col}_P"] = p_value
                results_df.loc[idx, f"{col}_D"] = cohen_d

    return results_df, common_indices, cc_columns


# def is_all_nan(lst):
#     if lst is np.nan:
#         return True
#     if isinstance(lst, list):
#         return all(np.isnan(x) for x in lst)
#     return False


def stats_exec(learning_xyz, tp_data):
    conditions = [x for x in learning_xyz if x != "[KEEP]"]

    results = {}

    for condition in conditions:
        print(condition + "...")
        results[condition] = {}
        results[condition]["metrics"] = pd.DataFrame(
            index=learning_xyz[condition]["x_full_df"].index
        )

        print("transfer w...")
        winner_list = []
        score_list = []
        for roundn in learning_xyz[condition]["w_full_prob_df"]:
            winner_list.append(
                learning_xyz[condition]["w_full_prob_df"][roundn]["SVM_winner"]
            )
            score_list.append(
                learning_xyz[condition]["w_full_prob_df"][roundn]["SVM_prob"]
            )

        combined_winner = pd.concat(winner_list, axis=1)
        results[condition]["metrics"]["SVM_winner"] = combined_winner.apply(
            reduce_to_single_column, axis=1
        )

        combined_score = pd.concat(score_list, axis=1)
        results[condition]["metrics"]["SVM_prob"] = combined_score.mean(axis=1)
        results[condition]["metrics"].loc[
            results[condition]["metrics"]["SVM_winner"].isna(), "SVM_prob"
        ] = np.nan

        print("calculate TPA...")
        TPA_list = []
        for replicate in tp_data[condition]:
            TPA_list.append(tp_data[condition][replicate])
        combined_TPA = pd.concat(TPA_list, axis=1)
        results[condition]["metrics"]["TPA"] = combined_TPA.mean(axis=1)
        results[condition]["metrics"]["dTPA"] = combined_TPA.std(axis=1)

        print("calculate CA...")
        results[condition]["metrics"]["CA_relevant"] = "no"
        results[condition]["class_abundance"] = {}
        classnames = []
        for classname in [
            classname
            for classname in results[condition]["metrics"]["SVM_winner"]
            if not pd.isnull(classname)
        ]:
            results_class = results[condition]["metrics"][
                (results[condition]["metrics"]["SVM_winner"] == classname)
                & (results[condition]["metrics"]["SVM_prob"] >= 0.95)
                & (~results[condition]["metrics"]["TPA"].isnull())
            ]
            classnames.append(classname)
            results[condition]["metrics"].loc[
                results_class.index, "CA_relevant"
            ] = "yes"
            results[condition]["class_abundance"][classname] = {}
            results[condition]["class_abundance"][classname]["mean"] = np.mean(
                results_class["TPA"]
            )
            results[condition]["class_abundance"][classname]["std"] = np.std(
                results_class["TPA"]
            )
            results[condition]["class_abundance"][classname]["count"] = len(
                results_class
            )

        print("calculate CC...")
        for classname in classnames:
            CC_list = []
            for roundn in learning_xyz[condition]["z_full_df"]:
                CC_list.append(
                    learning_xyz[condition]["z_full_df"][roundn][classname]
                )
            combined_CC = pd.concat(CC_list, axis=0)
            # combined_CC = pd.concat(CC_list, axis = 1)

            results[condition]["metrics"]["CC_" + classname] = (
                combined_CC.median(axis=0)
            )
            Q1 = combined_CC.quantile(0.25)
            Q3 = combined_CC.quantile(0.75)
            IQR = Q3 - Q1

            results[condition]["metrics"]["dCC_" + classname] = IQR

            # results[condition]['metrics']['CC_'+classname] = combined_CC.mean(axis = 1)
            # results[condition]['metrics']['dCC_'+classname] = combined_CC.std(axis = 1)

        print("calculate pCC...")
        for classname in classnames:
            results[condition]["metrics"]["pCC_" + classname] = (
                results[condition]["metrics"]["CC_" + classname]
                / results[condition]["class_abundance"][classname]["mean"]
            )
            results[condition]["metrics"]["dpCC_" + classname] = (
                abs(
                    1
                    / results[condition]["class_abundance"][classname]["mean"]
                )
                * results[condition]["metrics"]["dCC_" + classname]
                + abs(
                    -results[condition]["metrics"]["CC_" + classname]
                    / (
                        results[condition]["class_abundance"][classname][
                            "mean"
                        ]
                    )
                    ** 2
                )
                * results[condition]["class_abundance"][classname]["std"]
            )

        print("calculate nCC...")
        nCC_col = [
            col
            for col in results[condition]["metrics"].columns
            if col.startswith("pCC")
        ]
        for col in nCC_col:
            nCC_col = col.replace("pCC_", "nCC_")
            results[condition]["metrics"][nCC_col] = results[condition][
                "metrics"
            ][col]
        nCC_cols = [
            col
            for col in results[condition]["metrics"].columns
            if col.startswith("nCC_")
        ]
        results[condition]["metrics"][nCC_cols] = results[condition][
            "metrics"
        ][nCC_cols].div(
            results[condition]["metrics"][nCC_cols].sum(axis=1), axis=0
        )
        for classname in classnames:
            results[condition]["metrics"]["CC_factor_" + classname] = (
                results[condition]["metrics"]["nCC_" + classname]
                / results[condition]["metrics"]["pCC_" + classname]
            )
            results[condition]["metrics"]["dnCC_" + classname] = (
                results[condition]["metrics"]["dpCC_" + classname]
                * results[condition]["metrics"]["CC_factor_" + classname]
            )

        print("calculate CPA...")
        for classname in classnames:
            results[condition]["metrics"]["CPA_" + classname] = (
                results[condition]["metrics"]["CC_" + classname]
                * results[condition]["metrics"]["TPA"]
            )
            results[condition]["metrics"]["dCPA_" + classname] = (
                abs(results[condition]["metrics"]["TPA"])
                * results[condition]["metrics"]["dCC_" + classname]
                + abs(results[condition]["metrics"]["CC_" + classname])
                * results[condition]["metrics"]["dTPA"]
            )

        print("calculate rCPA...")
        for classname in classnames:
            results[condition]["metrics"]["rCPA_" + classname] = (
                results[condition]["metrics"]["nCC_" + classname]
                * results[condition]["metrics"]["TPA"]
            )
            results[condition]["metrics"]["drCPA_" + classname] = (
                abs(results[condition]["metrics"]["TPA"])
                * results[condition]["metrics"]["dnCC_" + classname]
                + abs(results[condition]["metrics"]["nCC_" + classname])
                * results[condition]["metrics"]["dTPA"]
            )

        results[condition]["classnames"] = classnames
    return results


def comp_exec(learning_xyz, results):
    conditions = [x for x in learning_xyz if x != "[KEEP]"]

    for condition in conditions:
        results[condition]["metrics"] = results[condition]["metrics"][
            ~results[condition]["metrics"].index.duplicated(keep="first")
        ]

    class_lists = []
    for condition in conditions:
        class_lists.append(results[condition]["classnames"])
    classnames = list(set(class_lists[0]).intersection(*class_lists[1:]))

    combinations = []
    for con_1 in conditions:
        for con_2 in conditions:
            if con_1 != con_2:
                combinations.append((con_1, con_2))

    for comb in combinations:
        results[comb] = {}
        results[comb]["intersection_data"] = pd.merge(
            learning_xyz[comb[0]]["x_full_df"],
            learning_xyz[comb[1]]["x_full_df"],
            left_index=True,
            right_index=True,
            how="inner",
        )
        results[comb]["comparison"] = pd.DataFrame(
            index=results[comb]["intersection_data"].index
        )

        metrics_own = results[comb[0]]["metrics"]
        metrics_other = results[comb[1]]["metrics"]

        common_idx = metrics_own.index.intersection(metrics_other.index)

        metrics_own = metrics_own.loc[common_idx]
        metrics_other = metrics_other.loc[common_idx]

        for classname in classnames:
            results[comb]["comparison"]["RL_" + classname] = (
                metrics_other["nCC_" + classname]
                - metrics_own["nCC_" + classname]
            )
            results[comb]["comparison"]["dRL_" + classname] = (
                metrics_other["dnCC_" + classname]
                + metrics_own["dnCC_" + classname]
            )

        RLS_columns = [
            col
            for col in results[comb]["comparison"].columns
            if col.startswith("RL_")
        ]
        dRLS_columns = [
            col
            for col in results[comb]["comparison"].columns
            if col.startswith("dRL_")
        ]

        results[comb]["comparison"]["RLS"] = (
            results[comb]["comparison"][RLS_columns].abs().sum(axis=1)
        )
        results[comb]["comparison"]["dRLS"] = (
            results[comb]["comparison"][dRLS_columns].abs().sum(axis=1)
        )

        for classname in classnames:
            results[comb]["comparison"]["RC_" + classname] = (
                metrics_other["CC_" + classname]
                - metrics_own["CC_" + classname]
            )
            results[comb]["comparison"]["dRC_" + classname] = (
                metrics_other["dCC_" + classname]
                + metrics_own["dCC_" + classname]
            )

        RCS_columns = [
            col
            for col in results[comb]["comparison"].columns
            if col.startswith("RC_")
        ]
        dRCS_columns = [
            col
            for col in results[comb]["comparison"].columns
            if col.startswith("dRC_")
        ]

        results[comb]["comparison"]["RCS"] = (
            results[comb]["comparison"][RCS_columns].abs().sum(axis=1)
        )
        results[comb]["comparison"]["dRCS"] = (
            results[comb]["comparison"][dRCS_columns].abs().sum(axis=1)
        )

        for classname in classnames:
            results[comb]["comparison"]["CPD_" + classname] = (
                metrics_other["rCPA_" + classname]
                - metrics_own["rCPA_" + classname]
            )
            results[comb]["comparison"]["dCPD_" + classname] = (
                metrics_other["drCPA_" + classname]
                + metrics_own["drCPA_" + classname]
            )
    return results


# filtered_columns = [col for col in results[comb]['comparison'].columns if col.startswith('RL_')]

# # Sum the absolute values of the filtered columns
# results[comb]['comparison']['RLS'] = results[comb]['comparison'][filtered_columns].abs().sum(axis=1)
