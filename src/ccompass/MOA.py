"""Multiple organelle analysis."""

import logging

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import ttest_ind

from ._utils import PrefixFilter, get_mp_ctx
from .core import KEEP, ComparisonModel, ResultsModel, XYZ_Model

logger = logging.getLogger(__package__)


def most_frequent_or_nan(row):
    """Return the most frequent value in a row, or np.nan if there's a tie."""
    counts = row.value_counts()
    # If the row is empty, return np.nan
    if counts.empty:
        return np.nan

    # If there's only one unique value in the row, return that value
    if len(counts) == 1:
        return counts.idxmax()

    # If the two most frequent values occur the same number of times, return np.nan
    if counts.iloc[0] == counts.iloc[1]:
        return np.nan

    return counts.idxmax()


def compare_lists(list1, list2):
    """Compute the sum of absolute differences between two lists.

    Ignores NaNs.
    """
    # Function to compare two lists and handle NaNs
    return sum(
        abs(a - b)
        for a, b in zip(list1, list2, strict=True)
        if not pd.isna(a) and not pd.isna(b)
    )


def is_all_nan(list_):
    return (
        all(np.isnan(x) for x in list_)
        if isinstance(list_, list)
        else np.isnan(list_)
    )


def perform_mann_whitney_t_tests_per_cell(
    df1: pd.DataFrame, df2: pd.DataFrame, prefix: str
) -> pd.DataFrame:
    """Perform Mann-Whitney U tests, t-tests, and compute Cohen's D.

    Groups are compared for each cell in the given DataFrames'
    columns starting with `prefix`. Those columns are assumed to
    contain lists of values for each group.
    """
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

    return results_df


def calculate_common_indices(df1, df2):
    return df1.index.intersection(df2.index)


def impute_data(x: pd.Series, s: float = 1.8, w: float = 0.3) -> pd.Series:
    """Impute missing values in a Series column using a normal distribution."""

    # Exclude 0s, compute mean and stddev
    x = x.replace(-np.inf, np.nan)
    mean = np.nanmean(x)
    std = np.nanstd(x)

    mean_imp = mean - s * std
    sigma = std * w

    nan_mask = x.isna()
    x.loc[nan_mask] = np.random.normal(mean_imp, sigma, nan_mask.sum())
    return x


def stats_proteome(
    learning_xyz: dict[str, XYZ_Model],
    fract_data,
    fract_conditions,
    reliability: float,
):
    """Proteome prediction / statistics."""
    logger.info("Performing proteome prediction...")
    conditions = [x for x in fract_conditions if x != KEEP]
    results = {}

    for condition in conditions:
        subcons = [x for x in learning_xyz if x.startswith(condition + "_")]

        combined_index = pd.DataFrame(index=[]).index
        for subcon in subcons:
            combined_index = combined_index.union(
                fract_data["class"][subcon].index
            )

        result = results[condition] = ResultsModel()
        result.metrics = pd.DataFrame(index=combined_index)

        ## add marker:
        result.metrics["marker"] = np.nan
        for subcon in subcons:
            marker_df = learning_xyz[subcon].W_train_df[
                ~learning_xyz[subcon].W_train_df.index.duplicated(keep="first")
            ]
            result.metrics["marker"] = result.metrics["marker"].fillna(
                marker_df
            )

        ## add SVM results:
        result.SVM["winner_combined"] = pd.DataFrame(
            index=result.metrics.index
        )
        result.SVM["prob_combined"] = pd.DataFrame(index=result.metrics.index)
        for subcon in subcons:
            xyz = learning_xyz[subcon]
            logger.info(f"Processing {condition} / {subcon}...")
            xyz.w_full_combined = pd.DataFrame(index=xyz.x_full_df.index)
            xyz.w_full_prob_combined = pd.DataFrame(index=xyz.x_full_df.index)

            for round_id, round_results in xyz.round_results.items():
                w_full_prob_df = round_results.w_full_prob_df
                xyz.w_full_combined = pd.merge(
                    xyz.w_full_combined,
                    w_full_prob_df["SVM_winner"],
                    left_index=True,
                    right_index=True,
                    how="left",
                )
                xyz.w_full_combined = xyz.w_full_combined.loc[
                    ~xyz.w_full_combined.index.duplicated(keep="first")
                ]
                xyz.w_full_combined.rename(
                    columns={"SVM_winner": f"{round_id}_SVM_winner"},
                    inplace=True,
                )
                xyz.w_full_prob_combined = pd.merge(
                    xyz.w_full_prob_combined,
                    w_full_prob_df["SVM_prob"],
                    left_index=True,
                    right_index=True,
                    how="left",
                )
                xyz.w_full_prob_combined = xyz.w_full_prob_combined.loc[
                    ~xyz.w_full_prob_combined.index.duplicated(keep="first")
                ]
                xyz.w_full_prob_combined.rename(
                    columns={"SVM_prob": f"{round_id}_SVM_prob"}, inplace=True
                )

            SVM_equal = xyz.w_full_combined.apply(
                lambda row: row.nunique() == 1, axis=1
            )
            xyz.w_full_combined["SVM_winner"] = np.where(
                SVM_equal,
                xyz.w_full_combined.iloc[:, 0],
                np.nan,
            )
            xyz.w_full_prob_combined["SVM_prob"] = (
                xyz.w_full_prob_combined.mean(axis=1)
            )

            xyz.w_full_combined = pd.merge(
                xyz.w_full_combined,
                xyz.w_full_prob_combined[["SVM_prob"]],
                left_index=True,
                right_index=True,
                how="left",
            )
            xyz.w_full_combined = xyz.w_full_combined.loc[
                ~xyz.w_full_combined.index.duplicated(keep="first")
            ]
            xyz.w_full_combined.loc[
                xyz.w_full_combined["SVM_winner"].isna(),
                "SVM_prob",
            ] = np.nan

            result.SVM["winner_combined"] = pd.merge(
                result.SVM["winner_combined"],
                xyz.w_full_combined["SVM_winner"].rename(
                    f"SVM_winner_{subcon}"
                ),
                left_index=True,
                right_index=True,
                how="left",
            )
            result.SVM["winner_combined"] = result.SVM["winner_combined"].loc[
                ~result.SVM["winner_combined"].index.duplicated(keep="first")
            ]
            result.SVM["prob_combined"] = pd.merge(
                result.SVM["prob_combined"],
                xyz.w_full_combined["SVM_prob"].rename(f"SVM_prob_{subcon}"),
                left_index=True,
                right_index=True,
                how="left",
            )
            result.SVM["prob_combined"] = result.SVM["prob_combined"].loc[
                ~result.SVM["prob_combined"].index.duplicated(keep="first")
            ]

        SVM_equal = result.SVM["winner_combined"].apply(
            lambda row: row.nunique() == 1, axis=1
        )
        # SVM_equal =
        SVM_major = result.SVM["winner_combined"].apply(
            most_frequent_or_nan, axis=1
        )
        SVM_major.name = "SVM_subwinner"
        result.metrics["SVM_winner"] = np.where(
            SVM_equal,
            result.SVM["winner_combined"].iloc[:, 0],
            np.nan,
        )
        result.metrics = pd.merge(
            result.metrics,
            SVM_major,
            left_index=True,
            right_index=True,
            how="left",
        )
        prob_means = result.SVM["prob_combined"].mean(axis=1)
        result.metrics["SVM_prob"] = np.nan
        result.metrics.loc[
            result.metrics["SVM_winner"].notna(), "SVM_prob"
        ] = prob_means

        ## add CClist:
        for subcon in subcons:
            #: TODO(performance): unnecessary conversions
            stacked_arrays = np.stack(
                list(
                    sr.z_full
                    for rr in learning_xyz[subcon].round_results.values()
                    for sr in rr.subround_results.values()
                )
            )
            learning_xyz[subcon].z_full_mean_df = pd.DataFrame(
                np.mean(stacked_arrays, axis=0),
                index=learning_xyz[subcon].x_full_df.index,
                columns=learning_xyz[subcon].Z_train_df.columns,
            )

        result.classnames = list(
            set(
                classname
                for subcon in subcons
                for classname in learning_xyz[subcon].classes
            )
        )

        for classname in result.classnames:
            CC_list = pd.DataFrame(index=combined_index)
            for subcon in subcons:
                CC_list = pd.merge(
                    CC_list,
                    learning_xyz[subcon]
                    .z_full_mean_df[classname]
                    .rename(f"{classname}_{subcon}"),
                    left_index=True,
                    right_index=True,
                    how="left",
                )

            CC_list["CClist_" + classname] = CC_list.apply(
                lambda row: row.tolist(), axis=1
            )

            result.metrics = pd.merge(
                result.metrics,
                CC_list["CClist_" + classname],
                left_index=True,
                right_index=True,
                how="left",
            )
            result.metrics = result.metrics.loc[
                ~result.metrics.index.duplicated(keep="first")
            ]

        ## add CC:
        for classname in result.classnames:
            result.metrics["CC_" + classname] = result.metrics[
                "CClist_" + classname
            ].apply(lambda x: np.nanmean(x) if x else np.nan)
        cc_cols = [
            col for col in result.metrics.columns if col.startswith("CC_")
        ]
        cc_sums = result.metrics[cc_cols].sum(axis=1, skipna=True)
        result.metrics[cc_cols] = result.metrics[cc_cols].div(cc_sums, axis=0)

        ## add NN_winner:
        cc_columns = result.metrics[
            [col for col in result.metrics.columns if col.startswith("CC_")]
        ]
        max_col = cc_columns.idxmax(axis=1)
        max_col = max_col.astype(str)
        result.metrics["NN_winner"] = max_col.str.replace("CC_", "")

        # add fCC:
        for class_act in result.classnames:
            nonmarker_z = result.metrics.loc[
                (result.metrics["marker"] != class_act)
                & (result.metrics["marker"].isna() == False)
            ][["CC_" + class_act]]
            thresh = np.percentile(
                nonmarker_z["CC_" + class_act].tolist(),
                reliability,
            )
            result.metrics["fCC_" + class_act] = result.metrics[
                "CC_" + class_act
            ]
            result.metrics.loc[
                result.metrics["fCC_" + class_act] < thresh,
                "fCC_" + class_act,
            ] = 0.0

        fcc_cols = [col for col in result.metrics if col.startswith("fCC_")]
        fcc_sums = result.metrics[fcc_cols].sum(axis=1)
        result.metrics[fcc_cols] = result.metrics[fcc_cols].div(
            fcc_sums, axis=0
        )

        ## add fNN_winner:
        fcc_columns = result.metrics[
            [col for col in result.metrics.columns if col.startswith("fCC_")]
        ]
        # compute idxmax(axis=1), but only for rows that aren't all NaN
        #  (The behavior of DataFrame.idxmax with all-NA values,
        #   or any-NA and skipna=False, is deprecated.)
        fmax_col = pd.Series(
            "nan",
            index=fcc_columns.index,
        )
        not_all_na_mask = ~fcc_columns.isna().all(axis=1)
        fmax_col[not_all_na_mask] = fcc_columns[not_all_na_mask].idxmax(axis=1)
        result.metrics["fNN_winner"] = fmax_col.str.replace("fCC_", "")

    logger.info("Proteome prediction done.")

    return results


def global_comparisons(
    results: dict[str, ResultsModel],
    max_processes: int = 1,
) -> dict[tuple[str, str], ComparisonModel]:
    """Compute global changes."""
    logger.info("Calculating global changes...")

    conditions = list(results)

    # deduplicate indices
    for result in results.values():
        result.metrics = result.metrics[
            ~result.metrics.index.duplicated(keep="first")
        ]

    arg_lists = [
        (con_1, con_2, results[con_1], results[con_2])
        for con_1 in conditions
        for con_2 in conditions
        if con_1 != con_2
    ]

    ctx = get_mp_ctx()
    with ctx.Pool(processes=max_processes) as pool:
        comparisons = dict(pool.map(_global_comparison_entry, arg_lists))

    logger.info("Global changes calculated.")

    return comparisons


def _global_comparison_entry(args):
    cond1, cond2, result1, result2 = args
    log_prefix = f"[{cond1} vs. {cond2}]"
    sub_logger = logger.getChild(log_prefix)
    sub_logger.addFilter(PrefixFilter(log_prefix))
    sub_logger.info("Starting global comparison...")
    return (cond1, cond2), global_comparison(result1, result2, sub_logger)


def global_comparison(
    result1: ResultsModel,
    result2: ResultsModel,
    logger: logging.Logger = logger,
) -> ComparisonModel:
    """Perform a single global comparison."""
    classnames = list(set(result1.classnames) & set(result2.classnames))
    comparison = ComparisonModel()

    metrics_own = result1.metrics
    metrics_other = result2.metrics

    ## prepare data:
    comparison.intersection_data = pd.merge(
        metrics_own,
        metrics_other,
        left_index=True,
        right_index=True,
        how="inner",
    )
    comparison.metrics = pd.DataFrame(index=comparison.intersection_data.index)

    logger.info("performing t-tests...")

    for classname in classnames:
        comparison.metrics["RL_" + classname] = (
            metrics_other["CC_" + classname] - metrics_own["CC_" + classname]
        )

    rl_cols = [
        col for col in comparison.metrics.columns if col.startswith("RL_")
    ]
    comparison.metrics["RLS"] = comparison.metrics[rl_cols].abs().sum(axis=1)

    for classname in classnames:
        comparison.metrics["fRL_" + classname] = (
            metrics_other["fCC_" + classname] - metrics_own["fCC_" + classname]
        )
    frl_cols = [
        col for col in comparison.metrics.columns if col.startswith("fRL_")
    ]
    comparison.metrics["fRLS"] = comparison.metrics[frl_cols].abs().sum(axis=1)

    test_df = perform_mann_whitney_t_tests_per_cell(
        metrics_own, metrics_other, "CClist_"
    )
    common_indices = test_df.index
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
    comparison.metrics = pd.merge(
        comparison.metrics,
        test_df,
        left_index=True,
        right_index=True,
        how="left",
    )

    logger.info("calculate RLS lists...")
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

        cclists_own_transposed = [list(values) for values in zip(*cclists_own)]
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
    comparison.RLS_results = pd.Series(RLS_results)
    comparison.RLS_null = pd.Series(RLS_null)

    comparison.metrics["P(t)_RLS"] = np.nan
    comparison.metrics["P(u)_RLS"] = np.nan
    for index in comparison.metrics.index:
        if index in common_indices:
            # Perform the t-test
            stat, p_value = ttest_ind(
                comparison.RLS_results.loc[index],
                comparison.RLS_null.loc[index],
                nan_policy="omit",
            )
            comparison.metrics.loc[index, "P(t)_RLS"] = p_value
            if (
                is_all_nan(comparison.RLS_results.loc[index])
                or is_all_nan(comparison.RLS_null.loc[index])
                or len(set(comparison.RLS_results.loc[index])) == 1
                or len(set(comparison.RLS_null.loc[index])) == 1
            ):
                comparison.metrics.loc[index, "P(u)_RLS"] = pd.NA
            else:
                stat_u, p_value_u = stats.mannwhitneyu(
                    comparison.RLS_results.loc[index],
                    comparison.RLS_null.loc[index],
                    alternative="two-sided",
                )
                comparison.metrics.loc[index, "P(u)_RLS"] = p_value_u
        else:
            comparison.metrics.loc[index, "P(t)_RLS"] = pd.NA
            comparison.metrics.loc[index, "P(u)_RLS"] = pd.NA

    return comparison


def class_comparisons(
    tp_data: dict[str, pd.DataFrame],
    results: dict[str, ResultsModel],
    comparisons: dict[tuple[str, str], ComparisonModel],
) -> None:
    """Compute class-centric changes.

    :param tp_data: Total proteome data.
    :param results: Results from the multi-organelle analysis. Will be updated.
    :param comparisons: Global comparisons. Will be updated.
    """
    logger.info("Calculating class-centric changes...")

    for condition, result in results.items():
        compute_class_centric_changes(
            result=result, tp_data=tp_data[condition]
        )

    logger.info("comparing...")

    combinations = [
        (con1, con2) for con1 in results for con2 in results if con1 != con2
    ]

    ## create nRL and nRLS:
    for comb in combinations:
        class_centric_comparison(
            results[comb[0]],
            results[comb[1]],
            comparisons[comb],
        )

    logger.info("Class-centric changes calculated.")


def compute_class_centric_changes(
    result: ResultsModel, tp_data: pd.DataFrame
) -> None:
    """Compute class-centric changes."""
    ## add TPA:
    logger.info("creating total protein amount...")
    tp_nontrans = tp_data.map(lambda x: 2**x)
    TPA_list = [tp_nontrans[replicate] for replicate in tp_data]
    combined_TPA = pd.concat(TPA_list, axis=1)
    result.metrics["TPA"] = combined_TPA.mean(axis=1)
    result.metrics = result.metrics.loc[
        ~result.metrics.index.duplicated(keep="first")
    ]

    # add class abundance:
    logger.info("adding class abundance...")
    result.metrics["CA_relevant"] = "no"
    result.class_abundance = {}
    for classname in result.classnames:
        results_class = result.metrics[
            (result.metrics["NN_winner"] == classname)
            & (~result.metrics["TPA"].isnull())
        ]
        result.metrics.loc[results_class.index, "CA_relevant"] = "yes"
        result.class_abundance[classname] = {
            "CA": np.median(results_class["TPA"]),
            "count": len(results_class),
        }

    ## add nCClist:
    logger.info("adding nCClist...")
    for classname in result.classnames:
        result.metrics["nCClist_" + classname] = result.metrics[
            "CClist_" + classname
        ].apply(
            lambda lst: [
                x * result.class_abundance[classname]["CA"]
                if not np.isnan(x)
                else np.nan
                for x in lst
            ]
        )

    ## add normalized class contributions:
    logger.info("adding normalized class contributions...")
    for classname in result.classnames:
        result.metrics["nCC_" + classname] = (
            result.metrics["fCC_" + classname]
            * result.class_abundance[classname]["CA"]
        )
    # normalize:
    nCC_cols = [col for col in result.metrics if col.startswith("nCC_")]
    nCC_sums = result.metrics[nCC_cols].sum(axis=1)
    nCC_sums[nCC_sums == 0] = 1
    result.metrics[nCC_cols] = result.metrics[nCC_cols].div(nCC_sums, axis=0)

    ## add CPA, nCPA, ...
    logger.info("adding CPA...")
    for classname in result.classnames:
        cpa = result.metrics["CPA_" + classname] = (
            result.metrics["CC_" + classname] * result.metrics["TPA"]
        )

        ncpa = result.metrics["nCPA_" + classname] = (
            result.metrics["nCC_" + classname] * result.metrics["TPA"]
        )

        with np.errstate(divide="ignore", invalid="ignore"):
            log_cpa = result.metrics["CPA_log_" + classname] = np.log2(cpa)
            log_ncpa = result.metrics["nCPA_log_" + classname] = np.log2(ncpa)

        result.metrics[f"CPA_imp_{classname}"] = impute_data(log_cpa)
        result.metrics[f"nCPA_imp_{classname}"] = impute_data(log_ncpa)


def class_centric_comparison(
    result1: ResultsModel, result2: ResultsModel, comparison: ComparisonModel
) -> None:
    """Perform class-centric comparison based on the result for the two given
    conditions."""

    metrics_own = result1.metrics
    metrics_other = result2.metrics
    common_indices = calculate_common_indices(metrics_own, metrics_other)

    if set(result1.classnames) != set(result2.classnames):
        # below, the assumption is that the classnames are the same
        raise AssertionError("Classes do not match.")

    for classname in result1.classnames:
        comparison.metrics["nRL_" + classname] = (
            metrics_other["nCC_" + classname] - metrics_own["nCC_" + classname]
        )

    logger.info("calculating nRL values...")
    nrl_cols = [
        col for col in comparison.metrics.columns if col.startswith("nRL_")
    ]
    comparison.metrics["nRLS"] = comparison.metrics[nrl_cols].abs().sum(axis=1)

    nRLS_results = {}
    nRLS_null = {}
    for ID in common_indices:
        ncclists_own = [
            metrics_own.loc[ID, "nCClist_" + classname]
            for classname in result1.classnames
        ]
        ncclists_other = [
            metrics_other.loc[ID, "nCClist_" + classname]
            for classname in result2.classnames
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
    comparison.nRLS_results = pd.Series(nRLS_results)
    comparison.nRLS_null = pd.Series(nRLS_null)

    comparison.metrics["P(t)_nRLS"] = np.nan
    comparison.metrics["P(u)_nRLS"] = np.nan
    # TODO(performance): vectorize
    for index in comparison.metrics.index:
        if index in common_indices:
            # Perform the t-test
            stat, p_value = ttest_ind(
                comparison.nRLS_results.loc[index],
                comparison.nRLS_null.loc[index],
                nan_policy="omit",
            )
            comparison.metrics.loc[index, "P(t)_nRLS"] = p_value
            if (
                is_all_nan(comparison.nRLS_results.loc[index])
                or is_all_nan(comparison.nRLS_null.loc[index])
                or len(set(comparison.nRLS_results.loc[index])) == 1
                or len(set(comparison.nRLS_null.loc[index])) == 1
            ):
                comparison.metrics.loc[index, "P(u)_nRLS"] = pd.NA
            else:
                stat_u, p_value_u = stats.mannwhitneyu(
                    comparison.nRLS_results.loc[index],
                    comparison.nRLS_null.loc[index],
                    alternative="two-sided",
                )
                comparison.metrics.loc[index, "P(u)_nRLS"] = p_value_u
        else:
            comparison.metrics.loc[index, "P(t)_nRLS"] = pd.NA
            comparison.metrics.loc[index, "P(u)_nRLS"] = pd.NA

    logger.info("calculating CPA values...")

    for classname in result1.classnames:
        comparison.metrics["CFC_" + classname] = (
            metrics_other["CPA_imp_" + classname]
            - metrics_own["CPA_imp_" + classname]
        )

        comparison.metrics["nCFC_" + classname] = (
            metrics_other["nCPA_imp_" + classname]
            - metrics_own["nCPA_imp_" + classname]
        )
