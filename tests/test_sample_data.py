"""Tests related to the sample dataset
at https://zenodo.org/records/13901167."""

from pathlib import Path

import pytest

from ccompass.core import (
    IDENTIFIER,
    KEEP,
    NA,
    FractDataset,
    MarkerSet,
    SessionModel,
    TotalProtDataset,
    read_fract_table,
    read_marker_table,
    read_tp_table,
)
from ccompass.FDP import start_fract_data_processing

sample_data_dir = Path(__file__).parents[1] / "sample_data"

if not sample_data_dir.exists():
    pytest.skip("Sample data directory not found", allow_module_level=True)


def get_fract_data() -> FractDataset:
    data_file = sample_data_dir / "C-COMPASS_ProteomeInput.tsv"
    df = read_fract_table(data_file)
    assert df.shape == (8841, 156)
    infix = "_DH_SA_60minDIA1cvFAIMS_Regeneron_MouseLiver_Gradient_"
    suffix = ".htrms.PG.Quantity"
    table = [
        # col_id, condition, replicate, fraction
        ["PG.ProteinGroups", IDENTIFIER, NA, NA],
        ["PG.Genes", KEEP, NA, NA],
        [f"[11] 20220323{infix}04_Fr01{suffix}", "HFHF", 1, 1],
        [f"[12] 20220323{infix}04_Fr02{suffix}", "HFHF", 1, 2],
        [f"[13] 20220323{infix}04_Fr03{suffix}", "HFHF", 1, 3],
        [f"[14] 20220323{infix}04_Fr04{suffix}", "HFHF", 1, 4],
        [f"[15] 20220323{infix}04_Fr05{suffix}", "HFHF", 1, 5],
        [f"[16] 20220323{infix}04_Fr06{suffix}", "HFHF", 1, 6],
        [f"[17] 20220323{infix}04_Fr07{suffix}", "HFHF", 1, 7],
        [f"[18] 20220323{infix}04_Fr08{suffix}", "HFHF", 1, 8],
        [f"[19] 20220323{infix}04_Fr09{suffix}", "HFHF", 1, 9],
        [f"[20] 20220323{infix}04_Fr10{suffix}", "HFHF", 1, 10],
        [f"[21] 20220323{infix}04_Fr11{suffix}", "HFHF", 1, 11],
        [f"[22] 20220323{infix}04_Fr12{suffix}", "HFHF", 1, 12],
        [f"[23] 20220323{infix}04_Fr13{suffix}", "HFHF", 1, 13],
        [f"[24] 20220323{infix}04_Fr14{suffix}", "HFHF", 1, 14],
        [f"[25] 20220323{infix}05_Fr01{suffix}", "HFHF", 2, 1],
        [f"[26] 20220323{infix}05_Fr02{suffix}", "HFHF", 2, 2],
        [f"[27] 20220323{infix}05_Fr03{suffix}", "HFHF", 2, 3],
        [f"[28] 20220323{infix}05_Fr04{suffix}", "HFHF", 2, 4],
        [f"[29] 20220323{infix}05_Fr05{suffix}", "HFHF", 2, 5],
        [f"[30] 20220323{infix}05_Fr06{suffix}", "HFHF", 2, 6],
        [f"[31] 20220323{infix}05_Fr07{suffix}", "HFHF", 2, 7],
        [f"[32] 20220323{infix}05_Fr08{suffix}", "HFHF", 2, 8],
        [f"[33] 20220323{infix}05_Fr09{suffix}", "HFHF", 2, 9],
        [f"[34] 20220323{infix}05_Fr10{suffix}", "HFHF", 2, 10],
        [f"[35] 20220323{infix}05_Fr11{suffix}", "HFHF", 2, 11],
        [f"[36] 20220323{infix}05_Fr12{suffix}", "HFHF", 2, 12],
        [f"[37] 20220323{infix}05_Fr13{suffix}", "HFHF", 2, 13],
        [f"[38] 20220323{infix}05_Fr14{suffix}", "HFHF", 2, 14],
        [f"[39] 20220323{infix}06_Fr01{suffix}", "HFHF-fasted", 1, 1],
        [f"[40] 20220323{infix}06_Fr02{suffix}", "HFHF-fasted", 1, 2],
        [f"[41] 20220323{infix}06_Fr03{suffix}", "HFHF-fasted", 1, 3],
        [f"[42] 20220323{infix}06_Fr04{suffix}", "HFHF-fasted", 1, 4],
        [f"[43] 20220323{infix}06_Fr05{suffix}", "HFHF-fasted", 1, 5],
        [f"[44] 20220323{infix}06_Fr06{suffix}", "HFHF-fasted", 1, 6],
        [f"[45] 20220323{infix}06_Fr07{suffix}", "HFHF-fasted", 1, 7],
        [f"[46] 20220323{infix}06_Fr08{suffix}", "HFHF-fasted", 1, 8],
        [f"[47] 20220323{infix}06_Fr09{suffix}", "HFHF-fasted", 1, 9],
        [f"[48] 20220323{infix}06_Fr10{suffix}", "HFHF-fasted", 1, 10],
        [f"[49] 20220323{infix}06_Fr11{suffix}", "HFHF-fasted", 1, 11],
        [f"[50] 20220323{infix}06_Fr12{suffix}", "HFHF-fasted", 1, 12],
        [f"[51] 20220323{infix}06_Fr13{suffix}", "HFHF-fasted", 1, 13],
        [f"[52] 20220323{infix}06_Fr14{suffix}", "HFHF-fasted", 1, 14],
        [f"[53] 20220323{infix}07_Fr01{suffix}", "HFHF", 3, 1],
        [f"[54] 20220323{infix}07_Fr02{suffix}", "HFHF", 3, 2],
        [f"[55] 20220323{infix}07_Fr03{suffix}", "HFHF", 3, 3],
        [f"[56] 20220323{infix}07_Fr04{suffix}", "HFHF", 3, 4],
        [f"[57] 20220323{infix}07_Fr05{suffix}", "HFHF", 3, 5],
        [f"[58] 20220323{infix}07_Fr06{suffix}", "HFHF", 3, 6],
        [f"[59] 20220323{infix}07_Fr07{suffix}", "HFHF", 3, 7],
        [f"[60] 20220323{infix}07_Fr08{suffix}", "HFHF", 3, 8],
        [f"[61] 20220323{infix}07_Fr09{suffix}", "HFHF", 3, 9],
        [f"[62] 20220323{infix}07_Fr10{suffix}", "HFHF", 3, 10],
        [f"[63] 20220323{infix}07_Fr11{suffix}", "HFHF", 3, 11],
        [f"[64] 20220323{infix}07_Fr12{suffix}", "HFHF", 3, 12],
        [f"[65] 20220323{infix}07_Fr13{suffix}", "HFHF", 3, 13],
        [f"[66] 20220323{infix}07_Fr14{suffix}", "HFHF", 3, 14],
        [f"[67] 20220323{infix}08_Fr01{suffix}", "HFHF", 4, 1],
        [f"[68] 20220323{infix}08_Fr02{suffix}", "HFHF", 4, 2],
        [f"[69] 20220323{infix}08_Fr03{suffix}", "HFHF", 4, 3],
        [f"[70] 20220323{infix}08_Fr04{suffix}", "HFHF", 4, 4],
        [f"[71] 20220323{infix}08_Fr05{suffix}", "HFHF", 4, 5],
        [f"[72] 20220323{infix}08_Fr06{suffix}", "HFHF", 4, 6],
        [f"[73] 20220323{infix}08_Fr07{suffix}", "HFHF", 4, 7],
        [f"[74] 20220323{infix}08_Fr08{suffix}", "HFHF", 4, 8],
        [f"[75] 20220323{infix}08_Fr09{suffix}", "HFHF", 4, 9],
        [f"[76] 20220323{infix}08_Fr10{suffix}", "HFHF", 4, 10],
        [f"[77] 20220323{infix}08_Fr11{suffix}", "HFHF", 4, 11],
        [f"[78] 20220323{infix}08_Fr12{suffix}", "HFHF", 4, 12],
        [f"[79] 20220323{infix}08_Fr13{suffix}", "HFHF", 4, 13],
        [f"[80] 20220323{infix}08_Fr14{suffix}", "HFHF", 4, 14],
        [f"[81] 20220323{infix}09_Fr01{suffix}", "HFHF-fasted", 2, 1],
        [f"[82] 20220323{infix}09_Fr02{suffix}", "HFHF-fasted", 2, 2],
        [f"[83] 20220323{infix}09_Fr03{suffix}", "HFHF-fasted", 2, 3],
        [f"[84] 20220323{infix}09_Fr04{suffix}", "HFHF-fasted", 2, 4],
        [f"[85] 20220323{infix}09_Fr05{suffix}", "HFHF-fasted", 2, 5],
        [f"[86] 20220323{infix}09_Fr06{suffix}", "HFHF-fasted", 2, 6],
        [f"[87] 20220323{infix}09_Fr07{suffix}", "HFHF-fasted", 2, 7],
        [f"[88] 20220323{infix}09_Fr08{suffix}", "HFHF-fasted", 2, 8],
        [f"[89] 20220323{infix}09_Fr09{suffix}", "HFHF-fasted", 2, 9],
        [f"[90] 20220323{infix}09_Fr10{suffix}", "HFHF-fasted", 2, 10],
        [f"[91] 20220323{infix}09_Fr11{suffix}", "HFHF-fasted", 2, 11],
        [f"[92] 20220323{infix}09_Fr12{suffix}", "HFHF-fasted", 2, 12],
        [f"[93] 20220323{infix}09_Fr13{suffix}", "HFHF-fasted", 2, 13],
        [f"[94] 20220323{infix}09_Fr14{suffix}", "HFHF-fasted", 2, 14],
        [f"[95] 20220323{infix}10_Fr01{suffix}", "HFHF-fasted", 3, 1],
        [f"[96] 20220323{infix}10_Fr02{suffix}", "HFHF-fasted", 3, 2],
        [f"[97] 20220323{infix}10_Fr03{suffix}", "HFHF-fasted", 3, 3],
        [f"[98] 20220323{infix}10_Fr04{suffix}", "HFHF-fasted", 3, 4],
        [f"[99] 20220323{infix}10_Fr05{suffix}", "HFHF-fasted", 3, 5],
        [f"[100] 20220323{infix}10_Fr06{suffix}", "HFHF-fasted", 3, 6],
        [f"[101] 20220323{infix}10_Fr07{suffix}", "HFHF-fasted", 3, 7],
        [f"[102] 20220323{infix}10_Fr08{suffix}", "HFHF-fasted", 3, 8],
        [f"[103] 20220323{infix}10_Fr09{suffix}", "HFHF-fasted", 3, 9],
        [f"[104] 20220323{infix}10_Fr10{suffix}", "HFHF-fasted", 3, 10],
        [f"[105] 20220323{infix}10_Fr11{suffix}", "HFHF-fasted", 3, 11],
        [f"[106] 20220323{infix}10_Fr12{suffix}", "HFHF-fasted", 3, 12],
        [f"[107] 20220323{infix}10_Fr13{suffix}", "HFHF-fasted", 3, 13],
        [f"[108] 20220323{infix}10_Fr14{suffix}", "HFHF-fasted", 3, 14],
        [f"[109] 20220326{infix}01_Fr01{suffix}", "chow", 1, 1],
        [f"[110] 20220326{infix}01_Fr02{suffix}", "chow", 1, 2],
        [f"[111] 20220326{infix}01_Fr03{suffix}", "chow", 1, 3],
        [f"[112] 20220326{infix}01_Fr04{suffix}", "chow", 1, 4],
        [f"[113] 20220326{infix}01_Fr05{suffix}", "chow", 1, 5],
        [f"[114] 20220326{infix}01_Fr06{suffix}", "chow", 1, 6],
        [f"[115] 20220326{infix}01_Fr07{suffix}", "chow", 1, 7],
        [f"[116] 20220326{infix}01_Fr08{suffix}", "chow", 1, 8],
        [f"[117] 20220326{infix}01_Fr09{suffix}", "chow", 1, 9],
        [f"[118] 20220326{infix}01_Fr10{suffix}", "chow", 1, 10],
        [f"[119] 20220326{infix}01_Fr11{suffix}", "chow", 1, 11],
        [f"[120] 20220326{infix}01_Fr12{suffix}", "chow", 1, 12],
        [f"[121] 20220326{infix}01_Fr13{suffix}", "chow", 1, 13],
        [f"[122] 20220326{infix}01_Fr14{suffix}", "chow", 1, 14],
        [f"[123] 20220326{infix}02_Fr01{suffix}", "chow", 2, 1],
        [f"[124] 20220326{infix}02_Fr02{suffix}", "chow", 2, 2],
        [f"[125] 20220326{infix}02_Fr03{suffix}", "chow", 2, 3],
        [f"[126] 20220326{infix}02_Fr04{suffix}", "chow", 2, 4],
        [f"[127] 20220326{infix}02_Fr05{suffix}", "chow", 2, 5],
        [f"[128] 20220326{infix}02_Fr06{suffix}", "chow", 2, 6],
        [f"[129] 20220326{infix}02_Fr07{suffix}", "chow", 2, 7],
        [f"[130] 20220326{infix}02_Fr08{suffix}", "chow", 2, 8],
        [f"[131] 20220326{infix}02_Fr09{suffix}", "chow", 2, 9],
        [f"[132] 20220326{infix}02_Fr10{suffix}", "chow", 2, 10],
        [f"[133] 20220326{infix}02_Fr11{suffix}", "chow", 2, 11],
        [f"[134] 20220326{infix}02_Fr12{suffix}", "chow", 2, 12],
        [f"[135] 20220326{infix}02_Fr13{suffix}", "chow", 2, 13],
        [f"[136] 20220326{infix}02_Fr14{suffix}", "chow", 2, 14],
        [f"[137] 20220326{infix}03_Fr01{suffix}", "chow", 3, 1],
        [f"[138] 20220326{infix}03_Fr02{suffix}", "chow", 3, 2],
        [f"[139] 20220326{infix}03_Fr03{suffix}", "chow", 3, 3],
        [f"[140] 20220326{infix}03_Fr04{suffix}", "chow", 3, 4],
        [f"[141] 20220326{infix}03_Fr05{suffix}", "chow", 3, 5],
        [f"[142] 20220326{infix}03_Fr06{suffix}", "chow", 3, 6],
        [f"[143] 20220326{infix}03_Fr07{suffix}", "chow", 3, 7],
        [f"[144] 20220326{infix}03_Fr08{suffix}", "chow", 3, 8],
        [f"[145] 20220326{infix}03_Fr09{suffix}", "chow", 3, 9],
        [f"[146] 20220326{infix}03_Fr10{suffix}", "chow", 3, 10],
        [f"[147] 20220326{infix}03_Fr11{suffix}", "chow", 3, 11],
        [f"[148] 20220326{infix}03_Fr12{suffix}", "chow", 3, 12],
        [f"[149] 20220326{infix}03_Fr13{suffix}", "chow", 3, 13],
        [f"[150] 20220326{infix}03_Fr14{suffix}", "chow", 3, 14],
    ]
    d = FractDataset(df=df, table=table)
    assert d.id_col == "PG.ProteinGroups"
    return d


def get_tp_data() -> TotalProtDataset:
    data_file = sample_data_dir / "C-COMPASS_ProteomeInput.tsv"
    df = read_tp_table(data_file)
    assert df.shape == (8841, 156)
    infix = (
        " 20220406_DH_SA_60minDIA1cvFAIMS_Regeneron_MouseLiver_TotalProteome_"
    )
    suffix = ".htrms.PG.Quantity"
    table = [
        # col_id, condition
        ["PG.ProteinGroups", IDENTIFIER],
        [f"[1]{infix}01{suffix}", "chow"],
        [f"[2]{infix}02{suffix}", "chow"],
        [f"[3]{infix}03{suffix}", "chow"],
        [f"[4]{infix}04{suffix}", "HFHF"],
        [f"[5]{infix}05{suffix}", "HFHF"],
        [f"[6]{infix}06{suffix}", "HFHF-fasted"],
        [f"[7]{infix}07{suffix}", "HFHF"],
        [f"[8]{infix}08{suffix}", "HFHF"],
        [f"[9]{infix}09{suffix}", "HFHF-fasted"],
        [f"[10]{infix}10{suffix}", "HFHF-fasted"],
    ]
    return TotalProtDataset(df=df, table=table)


def get_markerlist() -> MarkerSet:
    file = sample_data_dir / "C-COMPASS_MarkerList.txt"
    df = read_marker_table(file)
    assert df.shape == (2789, 2)
    return MarkerSet(
        df=df, identifier_col="Genename", class_col="MarkerCompartment"
    )


def get_sample_session_input() -> SessionModel:
    sess = SessionModel()
    sess.fract_input = {"C-COMPASS_ProteomeInput.tsv": get_fract_data()}
    sess.tp_input = {"C-COMPASS_ProteomeInput.tsv": get_tp_data()}
    sess.marker_sets = {"C-COMPASS_MarkerList.txt": get_markerlist()}
    sess.marker_fractkey = "PG.Genes"
    return sess


def test_save_sample_data_session():
    """Create and save a session with the sample dataset.

    For testing, to avoid the manual metadata input.
    """
    sess = get_sample_session_input()
    sess.to_numpy(Path(__file__).parent / "_sample_data_input.npy")


def test_fdp():
    sess = get_sample_session_input()
    (
        data_ways,
        std_ways,
        protein_info,
        fract_conditions,
    ) = start_fract_data_processing(
        sess.fract_input,
        sess.fract_preparams,
    )

    assert fract_conditions == ["[KEEP]", "HFHF", "HFHF-fasted", "chow"]
    assert (
        protein_info["PG.Genes"].reset_index()["PG.Genes"]
        == sess.fract_input["C-COMPASS_ProteomeInput.tsv"].df["PG.Genes"]
    ).all()
