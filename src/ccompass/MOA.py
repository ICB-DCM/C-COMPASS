"""Multiple organelle analysis."""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import ttest_ind

from .core import NeuralNetworkParametersModel


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


def compare_lists(list1, list2):
    # Function to compare two lists and handle NaNs
    return sum(
        abs(a - b)
        for a, b in zip(list1, list2)
        if not pd.isna(a) and not pd.isna(b)
    )


def is_all_nan(list_):
    return (
        all(np.isnan(x) for x in list_)
        if isinstance(list_, list)
        else np.isnan(list_)
    )


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


def calculate_common_indices(df1, df2):
    return df1.index.intersection(df2.index)


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


def stats_proteome(
    learning_xyz,
    NN_params: NeuralNetworkParametersModel,
    fract_data,
    fract_conditions,
):
    conditions = [x for x in fract_conditions if x != "[KEEP]"]
    results = {}

    for condition in conditions:
        subcons = [x for x in learning_xyz if x.startswith(condition + "_")]

        # ---------------
        ### if mode == 'deep':
        ###     combined_index = tp_data[condition].index
        # ---------------

        # elif mode == 'rough':
        #     empty_index = []
        #     df = pd.DataFrame(index = empty_index)
        #     combined_index = df.index

        empty_index = []
        df = pd.DataFrame(index=empty_index)
        combined_index = df.index

        for subcon in subcons:
            combined_index = combined_index.union(
                fract_data["class"][subcon].index
            )

        results[condition] = {}
        results[condition]["metrics"] = pd.DataFrame(index=combined_index)

        # ---------------
        ### TPA_list = []
        ### ## add TPA:
        ### if mode == 'deep':
        ###     TPA_list = []
        ###     tp_nontrans = tp_data[condition].map(lambda x: 2 ** x)
        ###     for replicate in tp_data[condition]:
        ###         TPA_list.append(tp_nontrans[replicate])
        ###     combined_TPA = pd.concat(TPA_list, axis = 1)
        ###     results[condition]['metrics']['TPA'] = combined_TPA.mean(axis = 1)
        ###     results[condition]['metrics'] = results[condition]['metrics'].loc[~results[condition]['metrics'].index.duplicated(keep='first')]
        ### elif mode == 'rough':
        ###     pass
        # ---------------

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
                learning_xyz[subcon]["w_full_combined"]["SVM_winner"].rename(
                    f"SVM_winner_{subcon}"
                ),
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
                learning_xyz[subcon]["w_full_combined"]["SVM_prob"].rename(
                    f"SVM_prob_{subcon}"
                ),
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
                    learning_xyz[subcon]["z_full_mean_df"][classname].rename(
                        f"{classname}_{subcon}"
                    ),
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

        ## add NN_winner:
        cc_columns = results[condition]["metrics"][
            [
                col
                for col in results[condition]["metrics"].columns
                if col.startswith("CC_")
            ]
        ]
        max_col = cc_columns.idxmax(axis=1)
        max_col = max_col.astype(str)
        results[condition]["metrics"]["NN_winner"] = max_col.str.replace(
            "CC_", ""
        )

        # ---------------
        ### ## add CA:
        ### if mode == 'deep':
        ###     results[condition]['metrics']['CA_relevant'] = 'no'
        ###     results[condition]['class_abundance'] = {}
        ###     for classname in classnames:
        ###     # for classname in [classname for classname in results[condition]['metrics']['SVM_winner'] if not pd.isnull(classname)]:
        ###         results_class = results[condition]['metrics'][(results[condition]['metrics']['NN_winner'] == classname) &
        ###                                                       (~results[condition]['metrics']['TPA'].isnull())]
        ###         results[condition]['metrics'].loc[results_class.index, 'CA_relevant'] = 'yes'
        ###         results[condition]['class_abundance'][classname] = {}
        ###         results[condition]['class_abundance'][classname]['CA'] = np.median(results_class['TPA'])
        ###         results[condition]['class_abundance'][classname]['count'] = len(results_class)
        ### elif mode == 'rough':
        ###     pass
        # ---------------

        # add fCC:
        for class_act in classnames:
            nonmarker_z = results[condition]["metrics"].loc[
                (results[condition]["metrics"]["marker"] != class_act)
                & (results[condition]["metrics"]["marker"].isna() == False)
            ][["CC_" + class_act]]
            thresh = np.percentile(
                nonmarker_z["CC_" + class_act].tolist(),
                NN_params.reliability,
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

        # ---------------
        ## add nCClist:
        ### if mode == 'deep':
        ###     for classname in classnames:
        ###         results[condition]['metrics']['nCClist_'+classname] = results[condition]['metrics']['CClist_'+classname].apply(lambda lst: [x * results[condition]['class_abundance'][classname]['CA'] if not np.isnan(x) else np.nan for x in lst])
        ###
        ###     #nCClist_df = results[condition]['metrics'][]
        ###
        ###
        ###     ## add nCC:
        ###     for classname in classnames:
        ###         results[condition]['metrics']['nCC_'+classname] = results[condition]['metrics']['fCC_'+classname] * results[condition]['class_abundance'][classname]['CA']
        ###     # normalize:
        ###     nCC_cols = [col for col in results[condition]['metrics'] if col.startswith('nCC_')]
        ###     nCC_sums = results[condition]['metrics'][nCC_cols].sum(axis=1)
        ###     nCC_sums[nCC_sums == 0] = 1
        ###     results[condition]['metrics'][nCC_cols] = results[condition]['metrics'][nCC_cols].div(nCC_sums, axis=0)
        ###
        # ## add CPA
        ###     for classname in classnames:
        ###         results[condition]['metrics']['CPA_'+classname] = results[condition]['metrics']['CC_'+classname] * results[condition]['metrics']['TPA']
        ###     for classname in classnames:
        ###         results[condition]['metrics']['nCPA_'+classname] = results[condition]['metrics']['nCC_'+classname] * results[condition]['metrics']['TPA']
        ### elif mode == 'rough':
        ###     pass
        # ---------------

        results[condition]["classnames"] = classnames

    return results


def global_comparison(results):
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

        # --------------
        # ## create RL, nRL, RLS, and nRLS:
        # for classname in classnames:
        #     comparison[comb]['metrics']['RL_'+classname] = results[comb[1]]['metrics']['CC_'+classname] - results[comb[0]]['metrics']['CC_'+classname]
        #     if mode == 'deep':
        #         comparison[comb]['metrics']['nRL_'+classname] = results[comb[1]]['metrics']['nCC_'+classname] - results[comb[0]]['metrics']['nCC_'+classname]
        #     elif mode == 'rough':
        #         pass
        # --------------

        for classname in classnames:
            comparison[comb]["metrics"]["RL_" + classname] = (
                results[comb[1]]["metrics"]["CC_" + classname]
                - results[comb[0]]["metrics"]["CC_" + classname]
            )

        rl_cols = [
            col
            for col in comparison[comb]["metrics"].columns
            if col.startswith("RL_")
        ]
        comparison[comb]["metrics"]["RLS"] = (
            comparison[comb]["metrics"][rl_cols].abs().sum(axis=1)
        )

        # --------------
        # if mode == 'deep':
        #     nrl_cols = [col for col in comparison[comb]['metrics'].columns if col.startswith('nRL_')]
        #     comparison[comb]['metrics']['nRLS'] = comparison[comb]['metrics'][nrl_cols].abs().sum(axis=1)
        # elif mode == 'rough':
        #     pass
        # --------------

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

        test_df, common_indices, ncc_columns = (
            perform_mann_whitney_t_tests_per_cell(
                metrics_own, metrics_other, "CClist_"
            )
        )
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

        # --------------
        # if mode == 'deep':
        #     nRLS_results = {}
        #     nRLS_null = {}
        #     for ID in common_indices:
        #         ncclists_own = [metrics_own.loc[ID, 'nCClist_' + classname] for classname in classnames]
        #         ncclists_other = [metrics_other.loc[ID, 'nCClist_' + classname] for classname in classnames]

        #         ncclists_own_transposed = [list(values) for values in zip(*ncclists_own)]
        #         ncclists_other_transposed = [list(values) for values in zip(*ncclists_other)]

        #         nRLS_results[ID] = []
        #         nRLS_null[ID] = []

        #         for i in range(len(ncclists_own_transposed)):
        #             for j in range(i+1, len(ncclists_own_transposed)):
        #                 null_result = compare_lists(ncclists_own_transposed[i], ncclists_own_transposed[j])
        #                 nRLS_null[ID].append(null_result)
        #         for i in range(len(ncclists_other_transposed)):
        #             for j in range(i+1, len(ncclists_other_transposed)):
        #                 null_result = compare_lists(ncclists_other_transposed[i], ncclists_other_transposed[j])
        #                 nRLS_null[ID].append(null_result)

        #         for own_list in ncclists_own_transposed:
        #             for other_list in ncclists_other_transposed:
        #                 comparison_result = compare_lists(own_list, other_list)
        #                 nRLS_results[ID].append(comparison_result)
        #     comparison[comb]['nRLS_results'] = pd.Series(nRLS_results)
        #     comparison[comb]['nRLS_null'] = pd.Series(nRLS_null)

        #     comparison[comb]['metrics']['P(t)_nRLS'] = np.nan
        #     comparison[comb]['metrics']['P(u)_nRLS'] = np.nan
        #     for index in comparison[comb]['metrics'].index:
        #         if index in common_indices:
        #             # Perform the t-test
        #             stat, p_value = ttest_ind(comparison[comb]['nRLS_results'].loc[index], comparison[comb]['nRLS_null'].loc[index], nan_policy='omit')
        #             comparison[comb]['metrics'].loc[index, 'P(t)_nRLS'] = p_value
        #             if is_all_nan(comparison[comb]['nRLS_results'].loc[index]) or is_all_nan(comparison[comb]['nRLS_null'].loc[index]) or len(set(comparison[comb]['nRLS_results'].loc[index])) == 1 or len(set(comparison[comb]['nRLS_null'].loc[index])) == 1:
        #                 comparison[comb]['metrics'].loc[index, 'P(u)_nRLS'] = pd.NA
        #             else:
        #                 stat_u, p_value_u = stats.mannwhitneyu(comparison[comb]['nRLS_results'].loc[index], comparison[comb]['nRLS_null'].loc[index], alternative='two-sided')
        #                 comparison[comb]['metrics'].loc[index, 'P(u)_nRLS'] = p_value_u
        #         else:
        #             comparison[comb]['metrics'].loc[index, 'P(t)_nRLS'] = pd.NA
        #             comparison[comb]['metrics'].loc[index, 'P(u)_nRLS'] = pd.NA
        # elif mode == 'rough':
        #     pass
        # --------------

        print("calculate CPAs...")

        # --------------
        # if mode == 'deep':
        #     for classname in classnames:
        #         metrics_own['CPA_log_'+classname] = np.log2(metrics_own['CPA_'+classname])
        #         metrics_own = impute_data(metrics_own, 'CPA_log_'+classname, 'CPA_imp_'+classname)

        #         metrics_other['CPA_log_'+classname] = np.log2(metrics_other['CPA_'+classname])
        #         metrics_other = impute_data(metrics_other, 'CPA_log_'+classname, 'CPA_imp_'+classname)

        #         comparison[comb]['metrics']['CFC_'+classname] = metrics_other['CPA_imp_'+classname] - metrics_own['CPA_imp_'+classname]

        #         metrics_own['nCPA_log_'+classname] = np.log2(metrics_own['nCPA_'+classname])
        #         metrics_own = impute_data(metrics_own, 'nCPA_log_'+classname, 'nCPA_imp_'+classname)

        #         metrics_other['nCPA_log_'+classname] = np.log2(metrics_other['nCPA_'+classname])
        #         metrics_other = impute_data(metrics_other, 'nCPA_log_'+classname, 'nCPA_imp_'+classname)

        #         comparison[comb]['metrics']['nCFC_'+classname] = metrics_other['nCPA_imp_'+classname] - metrics_own['nCPA_imp_'+classname]
        # elif mode == 'rough':
        #     pass
        # --------------

    return comparison


def class_comparison(tp_data, fract_conditions, results, comparison):
    conditions = [x for x in results]

    for condition in conditions:
        # combined_index = tp_data[condition].index
        classnames = results[condition]["classnames"]

        print("creating TPA...")
        ## add TPA:
        TPA_list = []
        TPA_list = []
        tp_nontrans = tp_data[condition].map(lambda x: 2**x)
        for replicate in tp_data[condition]:
            TPA_list.append(tp_nontrans[replicate])
        combined_TPA = pd.concat(TPA_list, axis=1)
        results[condition]["metrics"]["TPA"] = combined_TPA.mean(axis=1)
        results[condition]["metrics"] = results[condition]["metrics"].loc[
            ~results[condition]["metrics"].index.duplicated(keep="first")
        ]

        print("adding CA...")
        ## add CA:
        results[condition]["metrics"]["CA_relevant"] = "no"
        results[condition]["class_abundance"] = {}
        for classname in classnames:
            results_class = results[condition]["metrics"][
                (results[condition]["metrics"]["NN_winner"] == classname)
                & (~results[condition]["metrics"]["TPA"].isnull())
            ]
            results[condition]["metrics"].loc[
                results_class.index, "CA_relevant"
            ] = "yes"
            results[condition]["class_abundance"][classname] = {}
            results[condition]["class_abundance"][classname]["CA"] = np.median(
                results_class["TPA"]
            )
            results[condition]["class_abundance"][classname]["count"] = len(
                results_class
            )

        print("adding CC...")
        ## add nCClist:
        for classname in classnames:
            results[condition]["metrics"]["nCClist_" + classname] = results[
                condition
            ]["metrics"]["CClist_" + classname].apply(
                lambda lst: [
                    x * results[condition]["class_abundance"][classname]["CA"]
                    if not np.isnan(x)
                    else np.nan
                    for x in lst
                ]
            )

        print("adding nCC...")
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

        print("adding CPA...")
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

    print("comparing...")

    combinations = []
    for con_1 in conditions:
        for con_2 in conditions:
            if con_1 != con_2:
                combinations.append((con_1, con_2))

    ## create , nRL, , and nRLS:

    for comb in combinations:
        metrics_own = results[comb[0]]["metrics"]
        metrics_other = results[comb[1]]["metrics"]
        common_indices = calculate_common_indices(metrics_own, metrics_other)

        for classname in classnames:
            comparison[comb]["metrics"]["nRL_" + classname] = (
                results[comb[1]]["metrics"]["nCC_" + classname]
                - results[comb[0]]["metrics"]["nCC_" + classname]
            )

        print("calculating nRL values...")
        nrl_cols = [
            col
            for col in comparison[comb]["metrics"].columns
            if col.startswith("nRL_")
        ]
        comparison[comb]["metrics"]["nRLS"] = (
            comparison[comb]["metrics"][nrl_cols].abs().sum(axis=1)
        )

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
                        ncclists_own_transposed[i], ncclists_own_transposed[j]
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
                comparison[comb]["metrics"].loc[index, "P(t)_nRLS"] = p_value
                if (
                    is_all_nan(comparison[comb]["nRLS_results"].loc[index])
                    or is_all_nan(comparison[comb]["nRLS_null"].loc[index])
                    or len(set(comparison[comb]["nRLS_results"].loc[index]))
                    == 1
                    or len(set(comparison[comb]["nRLS_null"].loc[index])) == 1
                ):
                    comparison[comb]["metrics"].loc[index, "P(u)_nRLS"] = pd.NA
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

        print("calculating CPA values...")

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
                metrics_other, "CPA_log_" + classname, "CPA_imp_" + classname
            )

            comparison[comb]["metrics"]["CFC_" + classname] = (
                metrics_other["CPA_imp_" + classname]
                - metrics_own["CPA_imp_" + classname]
            )

            metrics_own["nCPA_log_" + classname] = np.log2(
                metrics_own["nCPA_" + classname]
            )
            metrics_own = impute_data(
                metrics_own, "nCPA_log_" + classname, "nCPA_imp_" + classname
            )

            metrics_other["nCPA_log_" + classname] = np.log2(
                metrics_other["nCPA_" + classname]
            )
            metrics_other = impute_data(
                metrics_other, "nCPA_log_" + classname, "nCPA_imp_" + classname
            )

            comparison[comb]["metrics"]["nCFC_" + classname] = (
                metrics_other["nCPA_imp_" + classname]
                - metrics_own["nCPA_imp_" + classname]
            )

    return comparison


def class_reset(results, comparison):
    for condition in results:
        results[condition]["class_abundance"] = {}
        results[condition]["metrics"].drop(
            ["TPA", "CA_relevant"], axis=1, inplace=True
        )

        classnames = results[condition]["classnames"]
        for classname in classnames:
            results[condition]["metrics"].drop(
                [
                    "nCClist_" + classname,
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
    for comb in comparison:
        comparison[comb]["metrics"].drop(
            ["nRLS", "P(t)_nRLS", "P(u)_nRLS"], axis=1, inplace=True
        )
        comparison[comb] = {
            k: v
            for k, v in comparison[comb].items()
            if k not in ["nRLS_results", "nRLS_null"]
        }
        classnames = results[comb[0]]["classnames"]
        for classname in classnames:
            comparison[comb]["metrics"].drop(
                ["nRL_" + classname, "CFC_" + classname, "nCFC_" + classname],
                axis=1,
                inplace=True,
            )
    return results, comparison
