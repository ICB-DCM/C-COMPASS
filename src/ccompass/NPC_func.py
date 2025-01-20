# FIXME: this module is currently unused
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler


def pre_process():
    print("pre-process run")
    return


def create_condition_dataset(
    data_init, tables_init, conditions, idents
):  # creates dictionary containing data sorted by condition
    dataset = {}
    data = {}
    identifier = []
    for path in tables_init:
        identifier = list(
            set(identifier + list(data_init[path][idents[path]]))
        )
    for condition in conditions:
        data_new = pd.DataFrame(index=identifier)
        for path in tables_init:
            for sample in tables_init[path]:
                data = pd.DataFrame()
                if sample[1] == condition:
                    samplename = sample[0]
                    data[samplename] = data_init[path][sample[0]]
                    data.set_index(data_init[path][idents[path]], inplace=True)
                    data_new = pd.merge(
                        data_new,
                        data,
                        right_index=True,
                        left_index=True,
                        how="outer",
                    )
                    if condition == "[KEEP]":
                        if samplename + "_x" in data_new.columns:
                            for element in list(data_new.index):
                                if pd.isnull(
                                    data_new[samplename + "_x"][element]
                                ):
                                    data_new[samplename + "_x"][element] = (
                                        data_new[samplename + "_y"][element]
                                    )
                                if pd.isnull(
                                    data_new[samplename + "_y"][element]
                                ):
                                    data_new[samplename + "_y"][element] = (
                                        data_new[samplename + "_x"][element]
                                    )
                        data_new = data_new.T.drop_duplicates().T
                        data_new.rename(
                            {samplename + "_x": samplename},
                            axis=1,
                            inplace=True,
                        )
                    else:
                        data_new = data_new.rename(
                            columns={
                                samplename: "Fr."
                                + str(sample[3])
                                + "_"
                                + samplename
                                + "_Rep."
                                + str(sample[2])
                            }
                        )
        dataset[condition] = data_new
    return dataset


def create_params(values, data):  # summarizes parameters
    params = {}
    cons = []
    if values["--scale1_condition--"]:
        params["pre-scaling"] = "conditions"
    elif values["--scale1_replicate--"]:
        params["pre-scaling"] = "replicates"
    else:
        params["pre-scaling"] = False
    params["min VV"] = values["--missing_number--"]
    params["min count"] = values["--min_count--"]
    if values["--correl_rep_active--"]:
        params["excluded worst"] = True
    else:
        params["excluded worst"] = False
    if values["--median--"]:
        params["mode"] = "median"
    elif values["--concat--"]:
        params["mode"] = "concatenated"
    elif values["--separate--"]:
        params["mode"] = "separate"
    if values["--scale2_0--"]:
        params["post-scaling"] = "MinMax (0-1)"
    elif values["--scale2_1--"]:
        params["post-scaling"] = "A=1"
    else:
        params["post-scaling"] = False
    if values["--zeros_active--"]:
        params["baselines removed"] = True
    else:
        params["baselines removed"] = False
    for condition in data:
        cons.append(condition)
    params["conditions"] = cons
    params["marker_identifier"] = "[IDENTIFIER]"
    return params


def scale_data(dataset):  # scales a dataset from 0 to 1
    scaler = MinMaxScaler()
    dataset_scaled = pd.DataFrame(
        scaler.fit_transform(dataset.T).T, columns=dataset.columns
    )
    dataset_scaled = dataset_scaled.set_index(dataset.index)
    return dataset_scaled


def create_repdata(
    condata,
):  # creates dictionary containing data sorted by replicate
    repdata = {}
    replist = []
    for sample in condata.columns:
        rep = sample[sample.rfind("_") + 1 :]
        if rep not in replist:
            replist.append(rep)
    for rep in replist:
        data = pd.DataFrame(index=condata.index)
        for sample in condata.columns:
            suffix = sample[sample.rfind("_") + 1 :]
            if suffix == rep:
                data = pd.merge(
                    data, condata[sample], left_index=True, right_index=True
                )
        repdata[rep] = data
    return repdata


def removeElements(
    lst, k
):  # removes elements from list with a count less than k
    counted = Counter(lst)
    return [el for el in lst if counted[el] >= k]


def calculate_innercorr(condata, fracts_corr, protlist_con, condition):
    icorr = pd.DataFrame(index=protlist_con[condition])
    for replicate in condata:
        repdata_own = condata[replicate]
        for fract in repdata_own.columns:
            prefix = fract[: fract.find("_")]
            if prefix not in fracts_corr[condition]:
                repdata_own = repdata_own.drop([fract], axis=1)
        correls = pd.DataFrame(index=protlist_con[condition])
        for rep in condata:
            if not rep == replicate:
                repdata_other = condata[rep]
                for fract in repdata_other.columns:
                    prefix = fract[: fract.find("_")]
                    if prefix not in fracts_corr[condition]:
                        repdata_other = repdata_other.drop([fract], axis=1)
                correls[rep] = np.nan
                for ID in protlist_con[condition]:
                    if ID in repdata_own.index and ID in repdata_other.index:
                        profile_own = repdata_own.loc[ID].tolist()
                        profile_other = repdata_other.loc[ID].tolist()
                        corr = pearsonr(profile_own, profile_other)[0]
                        correls[rep][ID] = corr
        correls[replicate] = correls.mean(axis=1)
        icorr = pd.merge(
            icorr,
            correls[replicate],
            left_index=True,
            right_index=True,
            how="outer",
        )
    return icorr


def def_corr_IDs(condata, protlist_con, condition, mincount):
    corr_IDs = []
    for ID in protlist_con[condition]:
        count = 0
        for replicate in condata:
            if ID in condata[replicate].index:
                count = count + 1
        if count > mincount:
            corr_IDs.append(ID)
    return corr_IDs


# -----------------------------------------------------------------------------
## Global functions:
# -----------------------------------------------------------------------------


def create_dataset(data_all, tables_all, identifier, conditions, values_NPC):
    for i in ["", "[IDENTIFIER]"]:  # create dataset to work with
        if i in conditions:
            conditions.remove(i)
    data_con_abs = create_condition_dataset(
        data_all, tables_all, conditions, identifier
    )
    data_keep = {}
    if "[KEEP]" in data_con_abs:
        data_keep = data_con_abs["[KEEP]"]
        del data_con_abs["[KEEP]"]
    params = create_params(values_NPC, data_con_abs)
    return data_con_abs, data_keep, params


def pre_scaling(data_con, data_con_abs, mode):
    if mode == "condition":  # pre-scaling over conditions
        data_con_scaled = {}
        for condition in data_con_abs:
            condata_scaled = scale_data(data_con_abs[condition])
            data_con_scaled[condition] = condata_scaled
            data_rep = create_repdata(condata_scaled)
            data_con[condition] = data_rep
    elif mode == "replicate":  # pre-scaling over replicates
        for condition in data_con_abs:
            data_rep = create_repdata(data_con_abs[condition])
            for replicate in data_rep:
                data = data_rep[replicate]
                data_scaled = scale_data(data)
                data_rep[replicate] = data_scaled
            data_con[condition] = data_rep
    else:  # no pre-scaling
        for condition in data_con_abs:
            data_rep = create_repdata(data_con_abs[condition])
            data_con[condition] = data_rep
    return data_con


def filter_missing(data_con, stats, values_NPC):
    data_con_missing = {}
    for condition in data_con:  # filter by Missing Values / Valid Values
        condata = data_con[condition]
        for replicate in condata:
            repdata = condata[replicate]
            count_before = len(repdata)
            repdata.dropna(
                thresh=values_NPC["--missing_number--"], inplace=True
            )
            repdata.replace(np.nan, 0.0, inplace=True)
            count_after = len(repdata)
            stats["filtered"]["by ValidValues"][
                condition + "_" + replicate
            ] = count_after - count_before
            condata[replicate] = repdata
        data_con_missing[condition] = condata
        # data_con[condition] = condata
    return data_con_missing


def filter_empty(data_con):
    for condition in data_con:  # filter empty frofiles
        condata = data_con[condition]
        for replicate in condata:
            repdata = condata[replicate]
            repdata = repdata[~(repdata == 0).all(axis=1)]
            condata[replicate] = repdata
        data_con[condition] = condata
    return data_con


def filter_count(data_con, mincount, stats):  # filter by count over replicates
    protlist_con = {}
    data_con_count = {}
    for condition in data_con:
        condata = data_con[condition]
        peplist = []
        for replicate in condata:
            peplist = peplist + list(condata[replicate].index)
        peplist = list(set(removeElements(peplist, mincount)))
        for replicate in condata:
            repdata = condata[replicate]
            count_before = len(repdata)
            for index in list(repdata.index):
                if index not in peplist:
                    repdata.drop(index, axis=0, inplace=True)
            count_after = len(repdata)
            stats["filtered"]["by count"][condition + "_" + replicate] = (
                count_after - count_before
            )
            condata[replicate] = repdata
        data_con_count[condition] = condata
        # data_con[condition] = condata
        protlist_con[condition] = peplist
    return data_con_count, protlist_con, stats


def list_samples(data_con):
    fracts_con = {}  # find and list all fractions in each condition
    fracts_count = {}
    fracts_corr = {}
    for condition in data_con:
        condata = data_con[condition]
        fracts_count[condition] = {}
        fractions = []
        for replicate in condata:
            repdata = condata[replicate]
            for sample in list(repdata.columns):
                prefix = sample[: sample.find("_")]
                fractnumber = int(prefix[3:])
                if fractnumber not in fractions:
                    fractions.append(fractnumber)
                    fracts_count[condition][fractnumber] = 1
                else:
                    fracts_count[condition][fractnumber] = (
                        fracts_count[condition][fractnumber] + 1
                    )
        fractions = sorted(fractions)
        fracts_con[condition] = fractions
        fracts_corr[condition] = [
            "Fr." + str(k)
            for k, v in fracts_count[condition].items()
            if v == max(fracts_count[condition].values())
        ]
    return fracts_con, fracts_count, fracts_corr


def calculate_icorr(data_con, fracts_corr, protlist_con):
    icorr = {}
    for condition in data_con:
        icorr_sub = calculate_innercorr(
            data_con[condition], fracts_corr, protlist_con, condition
        )
        # icorr[condition] = ray.get(icorr_sub)
        icorr[condition] = icorr_sub
        icorr[condition].fillna(0.0, inplace=True)
    return icorr


def remove_worst(data_con, protlist_con, mincount, stats, icorr):
    check_IDs = {}
    for condition in data_con:
        corr_IDs = def_corr_IDs(
            data_con[condition], protlist_con, condition, mincount
        )
        check_IDs[condition] = corr_IDs
    for condition in data_con:
        for replicate in data_con[condition]:
            stats["filtered"]["by InnerCorrelation"][
                condition + "_" + replicate
            ] = 0
    data_con_cleaned = {}
    for condition in data_con:
        # condata = copy.deepcopy(data_con[condition])
        condata = data_con[condition]
        correls = icorr[condition]
        for ID in check_IDs[condition]:
            minrep = correls.idxmin(axis=1)[ID]
            try:
                condata[minrep] = condata[minrep].drop(ID, axis=0)
            except Exception:
                print("correls: ")
                print(correls)
                print("minrep: " + str(minrep))
                print(correls.loc[[ID]])
                print("ID: " + str(ID))
            stats["filtered"]["by InnerCorrelation"][
                condition + "_" + minrep
            ] -= 1
        data_con_cleaned[condition] = condata
    return data_con_cleaned


def implement_icorr(data_keep, icorr):
    for condition in icorr:
        col_new = "InnerCorrelation_" + condition
        correls = pd.DataFrame()
        correls[col_new] = icorr[condition].mean(axis=1)
        data_keep = pd.merge(
            data_keep, correls, left_index=True, right_index=True
        )
    return data_keep


def create_median(data_con, fracts_con, scale):
    data_con_median = {}
    data_con_std = {}
    for condition in data_con:
        condata = data_con[condition]
        fractions = fracts_con[condition]
        con_vals = pd.DataFrame()
        con_std = pd.DataFrame()
        for fract in fractions:
            fract_vals = pd.DataFrame()
            fract_std = pd.DataFrame()
            prefix = "Fr." + str(fract)
            for replicate in condata:
                repdata = condata[replicate]
                for sample in repdata:
                    if sample[: sample.find("_")] == prefix:
                        fract_vals = pd.merge(
                            fract_vals,
                            repdata[sample],
                            left_index=True,
                            right_index=True,
                            how="outer",
                        )
                        fract_std = pd.merge(
                            fract_std,
                            repdata[sample],
                            left_index=True,
                            right_index=True,
                            how="outer",
                        )
            cols = [col for col in fract_vals.columns]
            fract_vals[condition + "_median_" + prefix] = fract_vals[
                cols
            ].median(axis=1)
            fract_std[condition + "_std_" + prefix] = fract_std[cols].std(
                axis=1
            )
            con_vals = pd.merge(
                con_vals,
                fract_vals[condition + "_median_" + prefix],
                left_index=True,
                right_index=True,
                how="outer",
            ).fillna(0.0)
            con_std = pd.merge(
                con_std,
                fract_std[condition + "_std_" + prefix],
                left_index=True,
                right_index=True,
                how="outer",
            ).fillna(0.0)
            data_con_std[condition] = con_std
            if scale:
                con_vals_scaled = scale_data(con_vals)
                data_con_median[condition] = con_vals_scaled
            else:
                data_con_median[condition] = con_vals
    return data_con_median, data_con_std


def create_concat(data_con, scale):
    data_con_concat = {}
    for condition in data_con:
        condata = data_con[condition]
        con_vals = pd.DataFrame()
        for replicate in condata:
            repdata = condata[replicate]
            renamedict = {}
            for sample in repdata:
                oldname = sample
                newname = (
                    condition
                    + sample[sample.rfind("_") :]
                    + "_"
                    + sample[: sample.find("_")]
                )
                renamedict[oldname] = newname
            con_vals = pd.merge(
                con_vals,
                repdata,
                left_index=True,
                right_index=True,
                how="outer",
            )
            con_vals.rename(renamedict, axis="columns", inplace=True)
        if scale:
            con_vals = scale_data(con_vals)
        data_con_concat[condition] = con_vals.fillna(0.0)
    return data_con_concat


def create_separate(data_con):
    data_con_sep = {}
    for condition in data_con:
        condata = data_con[condition]
        for replicate in condata:
            new_con = condition + "_" + replicate
            data_con_sep[new_con] = condata[replicate]
    return data_con_sep


def scale_area(data_con):
    data_con_area = {}
    for condition in data_con:
        condata = data_con[condition]
        data_con_area[condition] = condata.div(condata.sum(axis=1), axis=0)
    return data_con_area


def remove_zeros(data_con, stats, params):
    for condition in data_con:
        condata = data_con[condition].fillna(0)
        count_before = len(condata)
        condata = condata[~(condata == 0).all(axis=1)]
        count_after = len(condata)
        stats["filtered"]["by baseline profiles"][
            condition + "_" + params["mode"]
        ] = count_after - count_before
        data_con[condition] = condata
    return data_con


def calculate_ocorr(data_con, data_keep):
    for condition in data_con:
        for con in data_con:
            if not con == condition:
                col_new = "OuterCorrelation_" + condition + "_" + con
                data_keep[col_new] = np.nan
                data_own = data_con[condition].fillna(0.0)
                data_other = data_con[con].fillna(0.0)
                # data_own = copy.deepcopy(data_con[condition]).fillna(0.)
                # data_other = copy.deepcopy(data_con[con]).fillna(0.)
                fracts_own = []
                fracts_other = []
                for fract in data_own.columns:
                    fracts_own.append(fract[fract.rfind("_") + 1 :])
                for fract in data_other.columns:
                    fracts_other.append(fract[fract.rfind("_") + 1 :])
                fracts_both = [x for x in fracts_own if x in fracts_other]
                for fract in data_own.columns:
                    suffix = fract[fract.rfind("_") + 1 :]
                    if suffix not in fracts_both:
                        data_own = data_own.drop([fract], axis=1)
                for fract in data_other.columns:
                    suffix = fract[fract.rfind("_") + 1 :]
                    if suffix not in fracts_both:
                        data_other = data_other.drop([fract], axis=1)
                for ID in data_own.index:
                    if ID in data_other.index:
                        profile_own = data_own.loc[ID].tolist()
                        profile_other = data_other.loc[ID].tolist()
                        corr = pearsonr(profile_own, profile_other)
                        data_keep[col_new][ID] = corr[0]
    return data_keep
