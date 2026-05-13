#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 14:25:28 2025

@author: charlie
"""
import sys, os
sys.path.append("/usr/local/lib")
# import lhapdf as lha
import numpy as np
import pandas as pd

path = os.path.dirname(os.path.abspath(__file__)) + "/"
print(path)

# import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

#import scipy.stats as st
#from scipy.interpolate import make_interp_spline



plt.rcParams["font.family"] = 'Arial'
# plt.rcParams["font.fantasy"] = 'xkcd'
# plt.style.use("seaborn-v0_8-talk")
# plt.style.use("fivethirtyeight")
plt.rcParams["text.usetex"] = True
plt.rcParams["hatch.linewidth"] = 1.2
#plt.set_cmap("Set1")


#Mp = 0.9383

obs = sys.argv[1]
n = sys.argv[2]
m = sys.argv[3]
experiment = sys.argv[4]

if experiment == "CLIC":
    Vs = "380"
elif experiment == "FCC-ee":
    Vs = "92"

if obs not in ["dsig", "ALL"]:
    print("Observable not available...quitting")
    sys.exit(0)
    
if obs == "dsig":
    ylabel = r"$\langle d\sigma \vert " + n + r"," + m + r"  \rangle$"


elif obs == "ALL":
    ylabel = r"$\langle A_{LL} \vert" + n + r"," + m + r"  \rangle$"

#if n == "0" or n == "1":
obs_to_plot = "<" + obs + "|" + n + ";" + m +">"
#if n == "2":
#    obs_to_plot = "<" + obs + "|+2,+" + m +">"


fill_hatches = ["/", "\\", "//", "\\\\", "x"]
fill_colors = ["red", "orange", "forestgreen"]
# linestyles = ["solid", "dashed", "dotted", "dashdot", (0, (5, 10)), (0, (3, 5, 1, 5, 1, 5))]


z1_edges = [["0.1", "0.2"], ["0.2", "0.3"], ["0.3", "0.4"],\
            ["0.4", "0.5"], ["0.5", "0.7"], ["0.7", "0.9"]]

#Q2_vals = ["4", "10", "50", "100", "500", "1000"]

Q2_vals = ["4", "10", "100"]

df_list_Q2 = {}

for Q2 in Q2_vals:
    
    df_list_z1bins = []
    
    for z1v in z1_edges:
        
        if experiment == "CLIC":
            filename = "CLIC/Fixed_Q2_" + Q2 + "_Vs_" + Vs + "_thetac_0.02_z1_" + z1v[0] + "_" + z1v[1] + "_separated.txt"
        elif experiment == "FCC-ee":
            filename = "FCC-ee/Fixed_Q2_" + Q2 + "_Vs_" + Vs + "_thetac_0.03_z1_" + z1v[0] + "_" + z1v[1] + "_separated.txt"
        df = pd.read_csv(filename, engine = "python")

        
        df_to_plot = pd.DataFrame()
        df_to_plot["z2"] = df.z2
        df_to_plot["obsU"] = df[obs_to_plot+"U"]
        df_to_plot["obsL"] = df[obs_to_plot+"L"]
        df_list_z1bins.append(df_to_plot)
    
    df_list_Q2[Q2] = df_list_z1bins
        


if obs == "dsig":
    y_sqrts = -0.02
    y_Q20 = -0.03
    y_zrange = 0.041
    y_U = -.02
    y_U_label = -.0225
    y_L = -.03
    y_L_label = -.0325
elif obs == "ALL":
    y_sqrts = -3e-3
    y_Q20 = -4.2e-3
    y_zrange = 4.e-3    
    y_U = y_sqrts
    y_L = y_Q20
    y_U_label = y_sqrts - 2.5e-4
    y_L_label = y_Q20- 2.5e-4

       
fig, axs = plt.subplots(nrows=2, ncols=3, figsize = (15,10), sharex = True, sharey = True)
plt.subplots_adjust(wspace = 0.0, hspace=0.0)

fig.text(0.05, 0.5, ylabel, va='center', ha='center', rotation='vertical', fontsize=30)

axes = axs.flatten()


for i in range(len(axes)):
    ax = axes[i]
    ax.set_xlim(0.05, 0.85)
    
    if obs == "dsig":
        ax.set_ylim(-.035,.0499)
    elif obs == "ALL":
        ax.set_ylim(-.00499,.00499)


    ax.tick_params("both", direction = "in", top = True, right = True, labelsize = 24, length = 5, which = "major")
    ax.tick_params("both", direction = "in", top = True, right = True, labelsize = 24, length = 3, which = "minor")
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator()) 
    ax.set_axisbelow(False)
    
    if i > 2:   
        ax.set_xlabel(r"$z_2$", size = 30)
        
    if i == 0:
        ax.text(.1, y_sqrts, r"$\sqrt{s} = " + Vs + "$ GeV", size = 24)
        ax.text(.1, y_Q20, r"$Q^2_0 = 3$ GeV$^2$", size = 24)

    ax.axhline(0, 1e-3, 5.1, color = 'black', lw = 1, zorder = 0)

        # ax.yaxis.set_major_formatter(tck.FuncFormatter(lambda x, _: '{:.0e}'.format(x)))

    z1v = z1_edges[i]
    z1_str = z1v[0]  + r" $ < z_1 <$ " + z1v[1]
    # ax.set_yscale('log')
    ax.text(.1, y_zrange, z1_str, size = 24)
    
    for j in range(len(Q2_vals)):
        
        Q2 = Q2_vals[j]
        Q2_label = r"$Q^2 = $ " + Q2 + r" GeV$^2$" 
              
        if i == 0:
            ax.plot(df_list_Q2[Q2][i]["z2"], df_list_Q2[Q2][i]["obsU"], lw = 2, ls = "solid",\
                    color = fill_colors[j],\
                label = Q2_label)
            ax.plot(df_list_Q2[Q2][i]["z2"], df_list_Q2[Q2][i]["obsL"], lw = 2, ls = "dashed",\
                            color = fill_colors[j])
            
        if i > 0:
            ax.plot(df_list_Q2[Q2][i]["z2"], df_list_Q2[Q2][i]["obsU"], lw = 2, ls = "solid", color = fill_colors[j])
            ax.plot(df_list_Q2[Q2][i]["z2"], df_list_Q2[Q2][i]["obsL"], lw = 2, ls = "dashed", color = fill_colors[j])

    if i == 0:
        fig.legend(loc = "outside upper center", bbox_to_anchor=(.5, .99), ncol = 3, frameon = False,  prop = {'size': 26})
        
    if i == 2:
        ax.plot([0.1,0.25], [y_U for el in [0.1,0.25]], lw = 2, ls = "solid",\
                color = "black",  label = r"$U$")
        ax.plot([0.1,0.25], [y_L for el in [0.1,0.25]], lw = 2, ls = "dashed",\
                color = "black", label = r"$L$")
        ax.text(.3,y_U_label, r"$U$", size = 24)
        ax.text(.3,y_L_label, r"$L$", size = 24)
        
        
    ax.text(0.5, 0.5, 'Preliminary', transform=ax.transAxes,
            fontsize=42, color='gray', alpha=0.4,
            ha='center', va='center', rotation=30)
        
plt.show()
fig.savefig(path + obs + "_" + n + "_" + m + "_separated_" + experiment + ".pdf", dpi = 400, bbox_inches = "tight")



