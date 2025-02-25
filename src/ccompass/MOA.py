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


def perform_mann_whitney_t_tests(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    subcons1: list[str],
    subcons2: list[str],
    class_names: list[str],
) -> pd.DataFrame:
    """Perform Mann-Whitney U tests, t-tests, and compute Cohen's d.

    Compute statistics for each class and protein.
    The samples are expected in df1["CC_{class}_{subcon}"]
    and df2["CC_{class}_{subcon}"] for each class in class_names
    and each subcon in subcons1 and subcons2.
    """
    merged = pd.merge(
        df1[
            [
                f"CC_{cls}_{subcon}"
                for subcon in subcons1
                for cls in class_names
            ]
        ],
        df2[
            [
                f"CC_{cls}_{subcon}"
                for subcon in subcons2
                for cls in class_names
            ]
        ],
        left_index=True,
        right_index=True,
        how="inner",
        suffixes=("#1", "#2"),
    )

    def get_stats(row):
        # Columns for Mann-Whitney U results, t-test results, and Cohen's d
        res = {
            f"{prefix}_{cls}": np.nan
            for prefix in ["U", "T", "D", "P(U)", "P(T)"]
            for cls in class_names
        }
        for cls in class_names:
            cols1 = [f"CC_{cls}_{subcon}" for subcon in subcons1]
            group1 = row[cols1]
            if group1.nunique() == 1:
                continue
            cols2 = [f"CC_{cls}_{subcon}" for subcon in subcons2]
            group2 = row[cols2]
            if group2.nunique() == 1:
                continue

            # Perform Mann-Whitney U test
            u_stat, p_value_u = stats.mannwhitneyu(
                group1, group2, alternative="two-sided"
            )

            # Perform t-test
            t_stat, p_value_t = stats.ttest_ind(
                group1, group2, equal_var=False, nan_policy="omit"
            )

            # Calculating Cohen's d
            diff = [
                value_1 - value_2 for value_1 in group1 for value_2 in group2
            ]
            mean_diff = np.mean(diff)
            std_diff = np.std(diff, ddof=1)
            cohen_d = mean_diff / std_diff if std_diff != 0 else np.nan

            # Storing results
            res[f"U_{cls}"] = u_stat
            res[f"P(U)_{cls}"] = p_value_u
            res[f"D_{cls}"] = cohen_d
            res[f"T_{cls}"] = t_stat
            res[f"P(T)_{cls}"] = p_value_t
        return pd.Series(res)

    results_df = merged.apply(get_stats, axis=1)
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
        result.subcons = subcons
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
        combine_svm_replicate_results(result, learning_xyz, subcons)
        svm_equal = result.SVM["winner_combined"].apply(
            lambda row: row.nunique() == 1, axis=1
        )
        svm_major = (
            result.SVM["winner_combined"]
            .apply(most_frequent_or_nan, axis=1)
            .rename("SVM_subwinner")
        )

        result.metrics["SVM_winner"] = np.where(
            svm_equal,
            result.SVM["winner_combined"].iloc[:, 0],
            np.nan,
        )
        result.metrics = pd.merge(
            result.metrics,
            svm_major,
            left_index=True,
            right_index=True,
            how="left",
        )
        result.metrics["SVM_prob"] = np.nan
        result.metrics.loc[
            result.metrics["SVM_winner"].notna(), "SVM_prob"
        ] = result.SVM["prob_combined"].mean(axis=1)

        ## add CClist:
        for subcon in subcons:
            # average NN outputs from different rounds and subrounds
            stacked_arrays = np.stack(
                list(
                    sr.z_full_df.values
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
            sorted(
                set(
                    classname
                    for subcon in subcons
                    for classname in learning_xyz[subcon].classes
                )
            )
        )

        # collect the CC0 values from the different replicates
        for classname in result.classnames:
            cc0_df = pd.DataFrame(index=combined_index)
            for subcon in subcons:
                cc0_df = pd.merge(
                    cc0_df,
                    learning_xyz[subcon]
                    .z_full_mean_df[classname]
                    .rename(f"CC_{classname}_{subcon}"),
                    left_index=True,
                    right_index=True,
                    how="left",
                )
            cc0_df[f"CC_{classname}"] = cc0_df.mean(axis=1, skipna=True)

            result.metrics = pd.merge(
                result.metrics,
                cc0_df,
                left_index=True,
                right_index=True,
                how="left",
            )
            result.metrics = result.metrics.loc[
                ~result.metrics.index.duplicated(keep="first")
            ]

        # TODO: replace CC (currently the NN output) by CC0 or similar,
        #    and fCC by CC to match the naming in the paper

        # normalize CC values to sum to 1
        cc_cols = [f"CC_{classname}" for classname in result.classnames]
        cc_sums = result.metrics[cc_cols].sum(axis=1, skipna=True)
        result.metrics[cc_cols] = result.metrics[cc_cols].div(cc_sums, axis=0)

        ## add NN_winner:
        max_col = result.metrics[cc_cols].idxmax(axis=1).astype(str)
        result.metrics["NN_winner"] = max_col.str.replace("CC_", "")

        # add fCC: filter out false positives, renormalize
        for class_act in result.classnames:
            cc = result.metrics["CC_" + class_act]
            nonmarker_z = cc[
                # markers for other classes
                (result.metrics["marker"] != class_act)
                & (result.metrics["marker"].isna() == False)
            ]
            thresh = np.percentile(
                nonmarker_z,
                reliability,
            )
            result.metrics["fCC_" + class_act] = cc.mask(cc < thresh, 0.0)

        fcc_cols = [col for col in result.metrics if col.startswith("fCC_")]
        fcc_sums = result.metrics[fcc_cols].sum(axis=1)
        result.metrics[fcc_cols] = result.metrics[fcc_cols].div(
            fcc_sums, axis=0
        )

        ## add fNN_winner:
        fcc_columns = result.metrics[fcc_cols]
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


def combine_svm_round_results(xyz: XYZ_Model):
    """Combine SVM result from the different training rounds.

    Combine the SVM predictions and probabilities from the different training
    rounds of a single replicate into a single DataFrame and compute the
    combined winner and mean probability.
    """

    # combine predictions
    xyz.w_full_combined = pd.DataFrame(index=xyz.x_full_df.index)

    for round_id, round_results in xyz.round_results.items():
        xyz.w_full_combined = pd.merge(
            xyz.w_full_combined,
            round_results.w_full_prob_df["SVM_winner"],
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

    svm_equal = xyz.w_full_combined.apply(
        lambda row: row.nunique() == 1, axis=1
    )
    xyz.w_full_combined["SVM_winner"] = np.where(
        svm_equal,
        xyz.w_full_combined.iloc[:, 0],
        np.nan,
    )
    xyz.w_full_combined = xyz.w_full_combined.loc[
        ~xyz.w_full_combined.index.duplicated(keep="first")
    ]

    # combine probabilities
    xyz.w_full_prob_combined = pd.DataFrame(index=xyz.x_full_df.index)

    for round_id, round_results in xyz.round_results.items():
        xyz.w_full_prob_combined = pd.merge(
            xyz.w_full_prob_combined,
            round_results.w_full_prob_df["SVM_prob"],
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

    xyz.w_full_prob_combined["SVM_prob"] = xyz.w_full_prob_combined.mean(
        axis=1
    )
    xyz.w_full_combined = pd.merge(
        xyz.w_full_combined,
        xyz.w_full_prob_combined[["SVM_prob"]],
        left_index=True,
        right_index=True,
        how="left",
    )
    xyz.w_full_combined.loc[
        xyz.w_full_combined["SVM_winner"].isna(),
        "SVM_prob",
    ] = np.nan


def combine_svm_replicate_results(
    result: ResultsModel,
    learning_xyz: dict[str, XYZ_Model],
    subcons: list[str],
):
    """Combine SVM results from the different replicates.

    Combine the SVM predictions and probabilities from the different
    replicates of a single condition into a single DataFrame and compute
    the combined winner and mean probability.

    Populates `result.SVM`.
    """
    result.SVM["winner_combined"] = pd.DataFrame(index=result.metrics.index)
    result.SVM["prob_combined"] = pd.DataFrame(index=result.metrics.index)
    for subcon in subcons:
        xyz = learning_xyz[subcon]
        logger.info(f"Processing {subcon}...")

        combine_svm_round_results(xyz)

        result.SVM["winner_combined"] = pd.merge(
            result.SVM["winner_combined"],
            xyz.w_full_combined["SVM_winner"].rename(f"SVM_winner_{subcon}"),
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

    if max_processes > 1:
        ctx = get_mp_ctx()
        with ctx.Pool(processes=max_processes) as pool:
            comparisons = dict(pool.map(_global_comparison_entry, arg_lists))
    else:
        comparisons = dict(
            _global_comparison_entry(args) for args in arg_lists
        )
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
    classnames = sorted(
        list(set(result1.classnames) & set(result2.classnames))
    )
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

    # calculate relocation = ΔCC
    for classname in classnames:
        comparison.metrics["RL_" + classname] = (
            metrics_other["CC_" + classname] - metrics_own["CC_" + classname]
        )

    # calculate relocalization scores = Σ_compartment|ΔCC_compartment|
    rl_cols = [
        col for col in comparison.metrics.columns if col.startswith("RL_")
    ]
    comparison.metrics["RLS"] = comparison.metrics[rl_cols].abs().sum(axis=1)

    # relocalization scores after filtering out false positives
    for classname in classnames:
        comparison.metrics["fRL_" + classname] = (
            metrics_other["fCC_" + classname] - metrics_own["fCC_" + classname]
        )
    frl_cols = [
        col for col in comparison.metrics.columns if col.startswith("fRL_")
    ]
    comparison.metrics["fRLS"] = comparison.metrics[frl_cols].abs().sum(axis=1)

    # add statistics
    test_df = perform_mann_whitney_t_tests(
        metrics_own,
        metrics_other,
        result1.subcons,
        result2.subcons,
        class_names=classnames,
    )
    common_indices = test_df.index

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
        # [[CC_class1_rep1, CC_class1_rep2, ...], [CC_class2_rep1, ...], ...]
        cclists_own = [
            metrics_own.loc[
                ID, [f"CC_{classname}_{subcon}" for subcon in result1.subcons]
            ].tolist()
            for classname in classnames
        ]
        cclists_other = [
            metrics_other.loc[
                ID, [f"CC_{classname}_{subcon}" for subcon in result2.subcons]
            ].tolist()
            for classname in classnames
        ]

        # [[CC_class1_rep1, CC_class2_rep1, ...], [CC_class1_rep2, ...], ...]
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
    RLS_results = pd.Series(RLS_results)
    RLS_null = pd.Series(RLS_null)

    comparison.metrics["P(t)_RLS"] = np.nan
    comparison.metrics["P(u)_RLS"] = np.nan
    for index in comparison.metrics.index:
        if index in common_indices:
            # Perform the t-test
            stat, p_value = ttest_ind(
                RLS_results.loc[index],
                RLS_null.loc[index],
                nan_policy="omit",
            )
            comparison.metrics.loc[index, "P(t)_RLS"] = p_value
            if (
                is_all_nan(RLS_results.loc[index])
                or is_all_nan(RLS_null.loc[index])
                or len(set(RLS_results.loc[index])) == 1
                or len(set(RLS_null.loc[index])) == 1
            ):
                comparison.metrics.loc[index, "P(u)_RLS"] = pd.NA
            else:
                stat_u, p_value_u = stats.mannwhitneyu(
                    RLS_results.loc[index],
                    RLS_null.loc[index],
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
        ca = result.class_abundance[classname]["CA"]
        for subcon in result.subcons:
            result.metrics[f"nCC_{classname}_{subcon}"] = (
                result.metrics[f"CC_{classname}_{subcon}"] * ca
            )

    ## add normalized class contributions:
    logger.info("adding normalized class contributions...")
    for classname in result.classnames:
        result.metrics["nCC_" + classname] = (
            result.metrics["fCC_" + classname]
            * result.class_abundance[classname]["CA"]
        )
    # normalize:
    nCC_cols = [f"nCC_{classname}" for classname in result.classnames]
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
            metrics_own.loc[
                ID, [f"nCC_{classname}_{subcon}" for subcon in result1.subcons]
            ].tolist()
            for classname in result1.classnames
        ]
        ncclists_other = [
            metrics_other.loc[
                ID, [f"nCC_{classname}_{subcon}" for subcon in result2.subcons]
            ].tolist()
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
    nRLS_results = pd.Series(nRLS_results)
    nRLS_null = pd.Series(nRLS_null)

    comparison.metrics["P(t)_nRLS"] = np.nan
    comparison.metrics["P(u)_nRLS"] = np.nan
    # TODO(performance): vectorize
    for index in comparison.metrics.index:
        if index in common_indices:
            # Perform the t-test
            stat, p_value = ttest_ind(
                nRLS_results.loc[index],
                nRLS_null.loc[index],
                nan_policy="omit",
            )
            comparison.metrics.loc[index, "P(t)_nRLS"] = p_value
            if (
                is_all_nan(nRLS_results.loc[index])
                or is_all_nan(nRLS_null.loc[index])
                or len(set(nRLS_results.loc[index])) == 1
                or len(set(nRLS_null.loc[index])) == 1
            ):
                comparison.metrics.loc[index, "P(u)_nRLS"] = pd.NA
            else:
                stat_u, p_value_u = stats.mannwhitneyu(
                    nRLS_results.loc[index],
                    nRLS_null.loc[index],
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
