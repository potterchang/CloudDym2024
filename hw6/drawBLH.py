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
casename = "pbl_half_PU_2" #evergreen"
filename = f"EddyFlux_{casename}-000000-000720.nc"

t0 = 60
t1 = 481
xIdx= 0
loc_list = ["Pasture","Urban","Pasture+Urban"]
ds = xr.open_dataset(f"{DatasetDir}/{filename}")
time = ds['time']
zc = ds['zc']
eddy = ds['wanothano'][t0:t1+1, :, xIdx]

toffset = 300
tsteplength = 2
tminute = toffset + t0 * tsteplength
HHMM = f"{str(tminute//60).zfill(2)}:{str(tminute%60).zfill(2)}"
time_range = pd.date_range(start=f"2024-01-01 {HHMM}", periods=(t1-t0+1), freq="2T")


#%% Draw Variables

fig, ax  = plt.subplots(3,1, figsize=[8,9], dpi=300, sharex=True, sharey=True)

ax[0].set_ylim([0,1500])
ax[0].minorticks_on()
#ax[0].set_ylabel('z [m]')
ax[1].set_ylabel('z [m]')
ax[-1].set_xlabel('Local Time')
ax[-1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax[0].set_title(f"{casename}", loc="right", fontsize=14, weight="bold")
ax[0].set_title(r"$\overline{w'\theta'}$ [K m/s]", loc="center", fontsize=15)
for i in range(3):
    xIdx = i
    eddy = ds['wanothano'][t0:t1+1, :, i]
    ax[i].set_title(f'{loc_list[xIdx]} Avg',loc='left')

    levels = np.arange(-0.02, 0.08, 0.005)#np.arange(0, 1.75, 0.2)
    cmap = plt.cm.get_cmap('RdYlBu_r', len(levels)+1)
    norm = BoundaryNorm(levels, ncolors=len(levels)+1, clip=False, extend='both')
    PC0 = ax[i].pcolormesh(time_range, zc, eddy.T, cmap=cmap, norm=norm)
    plt.colorbar(PC0)

plt.savefig(f"EddyFlux-{casename}.png", bbox_inches='tight')

