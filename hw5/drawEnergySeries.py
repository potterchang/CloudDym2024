#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 20:13:34 2024

@author: yhc2080
"""
#%% Initialize
# Basic
import os
import sys
import time
import numpy as np
import pandas as pd
#import geopandas as gpd
# I/O Processing
import pickle
from scipy.interpolate import interpn, RegularGridInterpolator
import xarray as xr
import netCDF4 as nc
import multiprocessing
from multiprocessing import Pool
# Taiwan VVM
sys.path.append("/data/yhc2080/UTIL")
from TaiwanVVMLoader import TaiwanVVMTOPO, TaiwanVVMData
# Visualization
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D

#%% Load Data
DatasetDir = "../DATA"
casename = "pbl_half_PU_uarea_2" #evergreen"
filename = f"Energy_{casename}-000000-000720.nc"

t0 = 60
t1 = 481
xIdx= 0
loc_list = ["Pasture","Urban"]

def loadData(t0, t1, xIdx):
    ds = xr.open_dataset(f"{DatasetDir}/{filename}")
    time = ds['time']
    wth = ds['wth'][t0:t1+1,  xIdx] * 1004
    wqv = ds['wqv'][t0:t1+1,xIdx] * 2.5e6
    fdswtoa = ds['fdswtoa'][t0:t1+1, xIdx]
    
    toffset = 300
    tsteplength = 2
    tminute = toffset + t0 * tsteplength
    HHMM = f"{str(tminute//60).zfill(2)}:{str(tminute%60).zfill(2)}"
    time_range = pd.date_range(start=f"2024-01-01 {HHMM}", periods=(t1-t0+1), freq="2T")    

    return wth, wqv, fdswtoa, time_range

wth1, wqv1, fdswtoa1, time_range = loadData(t0, t1, xIdx=0)
wth2, wqv2, fdswtoa2, time_range = loadData(t0, t1, xIdx=1)

#%% Draw Energy

fig, ax  = plt.subplots(1,1, figsize=[6,4], dpi=300, sharex=True, sharey=True)

ax.set_ylim([-5,800])
ax.set_ylabel('[W m-2]')
ax.set_xlabel('Local Time')
ax.set_xlim([time_range[0], time_range[-1]])
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.minorticks_on()

ax.set_title(f"Energy Series", loc="left", fontsize=14, weight="bold")
ax.set_title(f"{casename}", loc="right", fontsize=14, weight="bold")

L1=ax.plot(time_range, wth1, color='limegreen', lw=2.5, ls='--',alpha=0.7)
L1=ax.plot(time_range, wth2, color='limegreen', lw=2.5, ls='-')
L2=ax.plot(time_range, wqv1, color='tab:blue', lw=1.8, ls='--',alpha=0.7)
L2=ax.plot(time_range, wqv2, color='tab:blue', lw=1.8, ls='-')
L3=ax.plot(time_range, fdswtoa1, color='tab:orange', lw=1.8, ls='--',alpha=0.7)
L3=ax.plot(time_range, fdswtoa2, color='tab:orange', lw=1.8, ls='-')

L3legend = Line2D([0], [0], linestyle='--', color='#666666',linewidth=1.8)
L4legend = Line2D([0], [0], linestyle='-', color='k',linewidth=1.8)
# Combine all legend elements
all_handles = [L1[0], L2[0], L3[0], L3legend, L4legend]
all_labels = ['wth','wqv','fdswtoa', 'Pasture Avg.', 'Urban Avg.']

# Add legend to plot
ax.legend(all_handles, all_labels, loc=1, ncol=1,
          fontsize=10,borderpad=0.2,handleheight=0.9,handlelength=1.5,handletextpad=0.3,
          labelspacing=0.2,columnspacing=1.0,framealpha=0.90)


plt.savefig(f"Energy-{casename}.png", bbox_inches='tight')