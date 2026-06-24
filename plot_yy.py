#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import sys

import matplotlib
from matplotlib.lines import Line2D
matplotlib.rc('text',usetex=True)
import pylab as py
from matplotlib import colors
import matplotlib.gridspec as gridspec

import matplotlib.pyplot as plt
import matplotlib.ticker as tck


plt.rcParams["font.family"] = 'Arial'
plt.rcParams["text.usetex"] = True
plt.rcParams["hatch.linewidth"] = 1.2

matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

plt.rcParams["ytick.minor.visible"] = True

from input_yy import np, _Vs, _thetac,\
_model,\
_kTmin, _kTmax,\
_etaqmin, _etaqmax,\
_etaqbmin, _etaqbmax

# print(_Vs)
if _Vs > 5000:
    Vs_plot_str = r"$\sqrt{s_{\rm NN}}$ = " + str(_Vs) + r" GeV"
    out_name = "LHC_PbPb"
else:
    Vs_plot_str = r"$\sqrt{s}$ = " + str(_Vs) + r" GeV"
    out_name = "FCC-ee"

if (_etaqmax + _etaqmin) == 0:
    eta_plot_str = r"$|\eta_q,\,\eta_{\bar q}| \leq $ " + str(_etaqmax)
else:
    eta_plot_str = str(_etaqmin) + r" $ \leq |\eta_q,\,\eta_{\bar q}| \leq $ " + str(_etaqmax)

# path = "/hpe-gr4/sanedda/Azimuthal-moments_yy/"

kine_str_dict = {"Vs" : "_Vs_" + str(_Vs),\
                "kT" : "kT_max_" + str(_kTmax) + "_kT_min_" + str(_kTmin),\
                "thetac" : "_thetac_" + str(_thetac),\
                "etaq" : "_etaq_" + str(_etaqmin) + "_" + str(_etaqmax),\
                "etaqb" : "_etaqb_" + str(_etaqbmin) + "_" + str(_etaqbmax)
                }
for key in kine_str_dict.keys():
    out_name += kine_str_dict[key]

z_ranges = [[0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 0.5], [0.5, 0.7], [0.7, 0.9]]
z_ranges_str = [["0.1", "0.2"], ["0.2", "0.3"], ["0.3", "0.4"], ["0.4", "0.5"], ["0.5", "0.7"], ["0.7", "0.9"]]

df_dict = {}
df_list = []
df_list_JAM = []

q_low = 0.02275
q_high = 0.97725


for zr in z_ranges_str:
    zstr = "_z1_" + zr[0] + "_" + zr[1]
    name = kine_str_dict["kT"] + kine_str_dict["Vs"] + kine_str_dict["thetac"] + zstr + kine_str_dict["etaq"] + kine_str_dict["etaqb"] + "_U_L.csv"
    if "JAM3D-2022" in _model:
        JAM_name = _model + "_" + name
    print(name)

    df_list.append(pd.read_csv(name))
    df_list_JAM.append(pd.read_csv(JAM_name))
    # df_dict[zstr.replace("_z1", "z1")] = pd.read_csv(name)


fig, axs = plt.subplots(ncols = 3, nrows = 2, figsize=(19, 14), sharey = "row")#, sharex = "col")

fig.subplots_adjust(hspace = 0, wspace = 0)

# fig.text(0.5, 0.5, "Preliminary", fontsize=100, color="gray",\
#     alpha=0.3, ha="center", va="center", rotation=30, fontweight="bold",
#     transform=fig.transFigure, zorder=0)

k = 0
for i in range(len(axs)):
    for j in range(len(axs[i])):
        ax = axs[i][j]
        # if i == 1:
            # ax.set_xlim(0, 1)
        ax.tick_params(which = "minor", direction = "in", top = True, right = True, labelsize = 18, length = 3)
        ax.tick_params(which = "major", direction = "in", top = True, right = True, labelsize = 18, length = 6)
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.set_axisbelow(False)
        ax.axhline(linewidth=.5, color="black")
        ax.set_xlim(0.05, 0.85)

        df = df_list[k]
        df_JAM = df_list_JAM[k]
        
        ax.plot(df.z2, df["dsig_U_median"], color = "blue", lw = 2)
        ax.plot(df.z2, df["dsig_L_median"], color = "red", lw = 2)

        ax.fill_between(df.z2, df["dsig_U_min"], df["dsig_U_max"], facecolor = "none", hatch = "\\", edgecolor = "blue", label = r"$U$ (our)")
        ax.fill_between(df.z2, df["dsig_L_min"], df["dsig_L_max"], facecolor = "none", hatch = "/", edgecolor = "red", label = r"$L$ (our)")
        
        ax.plot(df_JAM.z2, df_JAM["dsig_U_median"], color = "orange", lw = 2)
        ax.plot(df_JAM.z2, df_JAM["dsig_L_median"], color = "green", lw = 2)

        ax.fill_between(df_JAM.z2, df_JAM["dsig_U_min"], df_JAM["dsig_U_max"], facecolor = "none", hatch = "/", edgecolor = "orange", label = r"$U$ (" + _model + ")")
        ax.fill_between(df_JAM.z2, df_JAM["dsig_L_min"], df_JAM["dsig_L_max"], facecolor = "none", hatch = "\\", edgecolor = "green", label = r"$L$ (" + _model + ")")

        ax.set_ylim(-0.3299, 0.3299)
        ax.annotate(r"$z_1 \in [$" + z_ranges_str[k][0] + r"$\,:\,$" + z_ranges_str[k][1] + r"$]$", (0.05, 0.1), xycoords = 'axes fraction', size = 24)

        if i == 1:
            ax.set_xlabel(r"$z_2$", size=24)
        if i == 0 and j == 2:
            ax.legend(fontsize=20, frameon=False)

        if i == 0 and j == 0:
            ax.annotate(Vs_plot_str, (0.05, 0.875), xycoords = "axes fraction", size = 24)
            ax.annotate(eta_plot_str,(0.05, 0.35), xycoords = "axes fraction", size = 24)
            ax.annotate(r"$k_T \in [$" + str(_kTmin) + r"$\,:\,$" + str(_kTmax) + r"$]$ GeV", (0.05, 0.225), xycoords = "axes fraction", size = 24)

        # if j == 0:
            # ax.set_ylabel(r"$\langle d\sigma^{\rm unp} | \cos\phi_{12}\rangle$", size = 24)

        k += 1
        # ax.text(xl, .07, xlabels_xsec_dict["etaj"] + r' $\in [$'+ str(kine_lims_dict["etaj"][0]) + r"$\,:\,$" + str(kine_lims_dict["etaj"][1]) + r"$]$", size = 24 )

    # elif round(_Vs, 0) >= 105.0:
        # ax.text(xl, -.055, Vs_label, size = 24)
        # ax.text(xl, -.0675,xlabels_xsec_dict["etaj"] + r' $\in [$'+ str(kine_lims_dict["etaj"][0]) + r"$\,:\,$" + str(kine_lims_dict["etaj"][1]) + r"$]$", size = 24 )

    fig.text(0.07, 0.5, r"$\langle d\sigma^{\rm unp} | \cos\phi_{12}\rangle$", va='center', rotation='vertical', size = 24)
    plt.savefig(out_name + "_U_L.pdf", bbox_inches = "tight", dpi=400)

plt.show()


fig, axs = plt.subplots(ncols = 3, nrows = 2, figsize=(19, 14), sharey = "row")#, sharex = "col")

fig.subplots_adjust(hspace = 0, wspace = 0)

# fig.text(0.5, 0.5, "Preliminary", fontsize=100, color="gray",\
#     alpha=0.3, ha="center", va="center", rotation=30, fontweight="bold",
#     transform=fig.transFigure, zorder=0)

k = 0
for i in range(len(axs)):
    for j in range(len(axs[i])):
        ax = axs[i][j]
        # if i == 1:
            # ax.set_xlim(0, 1)
        ax.tick_params(which = "minor", direction = "in", top = True, right = True, labelsize = 18, length = 3)
        ax.tick_params(which = "major", direction = "in", top = True, right = True, labelsize = 18, length = 6)
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.set_axisbelow(False)
        ax.axhline(linewidth=.5, color="black")
        ax.set_xlim(0.05, 0.85)

        df = df_list[k]
        df_JAM = df_list_JAM[k]

        ax.plot(df.z2, df["ALL_U_median"], color = "blue", lw = 2)
        ax.plot(df.z2, df["ALL_L_median"], color = "red", lw = 2)

        ax.fill_between(df.z2, df["ALL_U_min"], df["ALL_U_max"], facecolor = "none", hatch = "\\", edgecolor = "blue", label = r"$U$ (our)")
        ax.fill_between(df.z2, df["ALL_L_min"], df["ALL_L_max"], facecolor = "none", hatch = "/", edgecolor = "red", label = r"$L$ (our)")    
        

        ax.plot(df_JAM.z2, df_JAM["ALL_U_median"], color = "orange", lw = 2)
        ax.plot(df_JAM.z2, df_JAM["ALL_L_median"], color = "green", lw = 2)

        ax.fill_between(df_JAM.z2, df_JAM["ALL_U_min"], df_JAM["ALL_U_max"], facecolor = "none", hatch = "/", edgecolor = "orange", label = r"$U$ (" + _model + ")")
        ax.fill_between(df_JAM.z2, df_JAM["ALL_L_min"], df_JAM["ALL_L_max"], facecolor = "none", hatch = "\\", edgecolor = "green", label = r"$L$ (" + _model + ")")


        ax.set_ylim(-0.3299, 0.3299)
        ax.annotate(r"$z_1 \in [$" + z_ranges_str[k][0] + r"$\,:\,$" + z_ranges_str[k][1] + r"$]$", (0.05, 0.1), xycoords = 'axes fraction', size = 24)

        if i == 1:
            ax.set_xlabel(r"$z_2$", size=24)
        if i == 0 and j == 2:
            ax.legend(fontsize=20, frameon=False)

        if i == 0 and j == 0:
            ax.annotate(Vs_plot_str, (0.05, 0.875), xycoords = "axes fraction", size = 24)
            ax.annotate(eta_plot_str,(0.05, 0.35), xycoords = "axes fraction", size = 24)
            ax.annotate(r"$k_T \in [$" + str(_kTmin) + r"$\,:\,$" + str(_kTmax) + r"$]$ GeV", (0.05, 0.225), xycoords = "axes fraction", size = 24)

        # if j == 0:
            # ax.set_ylabel(r"$\langle d\sigma^{\rm unp} | \cos\phi_{12}\rangle$", size = 24)

        k += 1
        # ax.text(xl, .07, xlabels_xsec_dict["etaj"] + r' $\in [$'+ str(kine_lims_dict["etaj"][0]) + r"$\,:\,$" + str(kine_lims_dict["etaj"][1]) + r"$]$", size = 24 )

    # elif round(_Vs, 0) >= 105.0:
        # ax.text(xl, -.055, Vs_label, size = 24)
        # ax.text(xl, -.0675,xlabels_xsec_dict["etaj"] + r' $\in [$'+ str(kine_lims_dict["etaj"][0]) + r"$\,:\,$" + str(kine_lims_dict["etaj"][1]) + r"$]$", size = 24 )

    fig.text(0.07, 0.5, r"$\langle A^{LL} | \cos\phi_{12}\rangle$", va='center', rotation='vertical', size = 24)
    plt.savefig(out_name + "_ALL_U_L.pdf", bbox_inches = "tight", dpi=400)

plt.show()
