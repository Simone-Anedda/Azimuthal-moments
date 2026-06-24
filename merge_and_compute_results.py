#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import sys
import glob

from input_yy import np, _Vs, _thetac,\
_model,\
_kTmin, _kTmax,\
_etaqmin, _etaqmax,\
_etaqbmin, _etaqbmax,\
_folder

print(_Vs)

GeV2tomb = 3.89379e-1

# path = "/hpe-gr4/sanedda/Azimuthal-moments_yy/" + _folder + "/"

path = ""

kine_str_dict = {"Vs" : "_Vs_" + str(_Vs),\
                "kT" : "kT_max_" + str(_kTmax) + "_kT_min_" + str(_kTmin),\
                "thetac" : "_thetac_" + str(_thetac),\
                "etaq" : "_etaq_" + str(_etaqmin) + "_" + str(_etaqmax),\
                "etaqb" : "_etaqb_" + str(_etaqbmin) + "_" + str(_etaqbmax)
                }

# z_ranges = [[0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 0.5], [0.5, 0.7], [0.7, 0.9]]
z_ranges_str = [["0.1", "0.2"], ["0.2", "0.3"], ["0.3", "0.4"], ["0.4", "0.5"], ["0.5", "0.7"], ["0.7", "0.9"]]

# df_dict = {}

cols = ["z2", "<dsig|c12>_U", "<dsig|c12>_L", "<ALL|c12>_U", "<ALL|c12>_L"]

cols_dict_mean = {"<dsig|c12>_U" : "dsig_U_mean", "<dsig|c12>_L" : "dsig_L_mean", "<ALL|c12>_U" : "ALL_U_mean", "<ALL|c12>_L" : "ALL_L_mean"}
cols_dict_median = {"<dsig|c12>_U" : "dsig_U_median", "<dsig|c12>_L" : "dsig_L_median", "<ALL|c12>_U" : "ALL_U_median", "<ALL|c12>_L" : "ALL_L_median"}
cols_dict_min = {"<dsig|c12>_U" : "dsig_U_min", "<dsig|c12>_L" : "dsig_L_min", "<ALL|c12>_U" : "ALL_U_min", "<ALL|c12>_L" : "ALL_L_min"}
cols_dict_max = {"<dsig|c12>_U" : "dsig_U_max", "<dsig|c12>_L" : "dsig_L_max", "<ALL|c12>_U" : "ALL_U_max", "<ALL|c12>_L" : "ALL_L_max"}

q_low = 0.02275
q_high = 0.97725


for zr in z_ranges_str:
    zstr = "_z1_" + zr[0] + "_" + zr[1]
    name = kine_str_dict["kT"] + kine_str_dict["Vs"] + kine_str_dict["thetac"] + zstr + kine_str_dict["etaq"] + kine_str_dict["etaqb"] + "_iset_"
    print(name)

    if not "JAM3D-2022" in _model:
        df_list = []
        for el in range(1, 2000, 200):
            df = pd.read_csv(path + name + str(el) + "_" + str(el + 199) + "_U_L.txt")
            df_list.append(df)

        df = pd.concat(df_list, ignore_index = True)

    else:
        if _model == "JAM3D-2022":
            set_str = str(1) + "_" + str(466)
        elif _model == "JAM3D-2022-nolat":
            set_str = str(1) + "_" + str(230)
        df = pd.read_csv(path + _model + "_" + name + set_str + "_U_L.txt")


    df_cols = df.columns

    if "denU" in df_cols and "denU" not in cols:
        df.denU *= GeV2tomb
        cols.append("denU")
        cols_dict_mean["denU"] = "denU_mean"
        cols_dict_median["denU"] = "denU_median"
        cols_dict_min["denU"] = "denU_min"
        cols_dict_max["denU"] = "denU_max"
    if "denL" in df_cols and "denL"  not in cols:
        df.denL *= GeV2tomb
        cols.append("denL")
        cols_dict_mean["denL"] = "denL_mean"
        cols_dict_median["denL"] = "denL_median"
        cols_dict_min["denL"] = "denL_min"
        cols_dict_max["denL"] = "denL_max"

    df = df[cols].groupby("z2")
    df_new_mean = df.mean().rename(columns = cols_dict_mean)
    df_new_median = df.median().rename(columns = cols_dict_median)
    df_new_min = df.quantile(q_low).rename(columns = cols_dict_min)
    df_new_max = df.quantile(q_high).rename(columns = cols_dict_max)

    df_full = pd.concat([df_new_mean, df_new_median, df_new_min, df_new_max], axis = 1)

    df_full["z2"] = df.mean().index

    cols_to_save = [col for col in df_full.columns if not "." in col]
#    print(cols_to_save)

    df_full[cols_to_save].to_csv(_model + "_" + name.replace("_iset_", "") + "_U_L.csv", index = False)
