#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 16:20:46 2024

@author: yhc2080
"""

#%% Initialize
# Basic
import os
import sys
import time
import numpy as np
import pandas as pd
# I/O Processing
import xarray as xr
import netCDF4 as nc
import multiprocessing
from multiprocessing import Pool
# Visualization
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
#%% Load Data
DatasetDir = "../DATA"
casename = "pbl_up" #evergreen"
filename = f"Turbulence_{casename}-000000-000720.nc"

t0 = 60
t1 = 481
xIdx= 0
loc_list = ["Pasture","Urban"]

def loadData(t0, t1, xIdx):
    ds = xr.open_dataset(f"{DatasetDir}/{filename}")
    time = ds['time']
    zc = ds['zc']
    tke = ds['tke'][t0:t1+1, :, xIdx]
    enstrophy = ds['enstrophy'][t0:t1+1, :, xIdx]
    th = ds['th'][t0:t1+1, :, xIdx]
    
    toffset = 300
    tsteplength = 2
    tminute = toffset + t0 * tsteplength
    HHMM = f"{str(tminute//60).zfill(2)}:{str(tminute%60).zfill(2)}"
    time_range = pd.date_range(start=f"2024-01-01 {HHMM}", periods=(t1-t0+1), freq="2T")
    
    # Calculate conventional BLH
    ## BLH1:  thbot + 0.5K
    thbot = th[:,1] 
    P05K = thbot + 0.5
    ind_posi = np.array((th - P05K) >= 0)
    adj = np.sum(ind_posi, axis=1) == 0
    ind_BLH1 = np.argmax(ind_posi, axis=1) - adj
    ind_BLH1 = np.ma.array(ind_BLH1, mask=ind_BLH1<0)
    BLH1 = np.take(zc, ind_BLH1)
    
    ## BLH2: max gradient
    dthdz = np.gradient(th, zc, axis=1)
    ind_BLH2 = np.argmax(dthdz, axis=1)
    BLH2 = np.take(zc, ind_BLH2)

    return tke, enstrophy, BLH1, BLH2, time_range, zc

tke0, enstrophy0, BLH10, BLH20, time_range, zc = loadData(t0, t1, xIdx=0)
tke1, enstrophy1, BLH11, BLH21, time_range, zc = loadData(t0, t1, xIdx=1)

#%% Draw BLH

fig, ax  = plt.subplots(1,1, figsize=[6,4], dpi=300, sharex=True, sharey=True)

ax.set_ylim([0,1500])
ax.set_ylabel('z [m]')
ax.set_xlabel('Local Time')
ax.set_xlim([time_range[0], time_range[-1]])
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.minorticks_on()

ax.set_title(f"PBL Depth", loc="left", fontsize=14, weight="bold")
ax.set_title(f"{casename}", loc="right", fontsize=14, weight="bold")

PC0 = ax.contour(time_range, zc, tke0.T, levels=[0.2], colors=["tab:red"], linewidths=2.0, linestyles='--',alpha=0.7)
PC0 = ax.contour(time_range, zc, tke1.T, levels=[0.2], colors=["tab:red"], linewidths=2.0, linestyles='-')
PC0_legend = Line2D([0], [0], linestyle='-', color='tab:red',linewidth=2.0)

PC1 = ax.contour(time_range, zc, enstrophy0.T*1e5, levels=[2], colors=["tab:blue"], linestyles='--',alpha=0.7)
PC1 = ax.contour(time_range, zc, enstrophy1.T*1e5, levels=[2], colors=["tab:blue"], linestyles='-')
PC1_legend = Line2D([0], [0], linestyle='-', color='tab:blue',linewidth=1.5)

L1=ax.plot(time_range, BLH10, color='limegreen', lw=2.5, ls='--',alpha=0.7)
L1=ax.plot(time_range, BLH11, color='limegreen', lw=2.5, ls='-')
L2=ax.plot(time_range, BLH20, color='tab:orange', lw=1.8, ls='--',alpha=0.7)
L2=ax.plot(time_range, BLH21, color='tab:orange', lw=1.8, ls='-')

L3legend = Line2D([0], [0], linestyle='--', color='#666666',linewidth=1.8)
L4legend = Line2D([0], [0], linestyle='-', color='k',linewidth=1.8)
# Combine all legend elements
all_handles = [PC0_legend, PC1_legend, L1[0], L2[0], L3legend, L4legend]
all_labels = ['TKE=0.2', 'Enstrophy=2e-5', r'$\theta_{sfc}+0.5K$', r'max $d\theta/dz$', 'Pasture Avg.', 'Urban Avg.']

# Add legend to plot
ax.legend(all_handles, all_labels, loc=2, ncol=1,
          fontsize=10,borderpad=0.2,handleheight=0.9,handlelength=1.5,handletextpad=0.3,
          labelspacing=0.2,columnspacing=1.0,framealpha=0.90)


plt.savefig(f"PBL-{casename}.png", bbox_inches='tight')
