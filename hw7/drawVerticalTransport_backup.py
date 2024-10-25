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
casename = "pbl_half_PU_1" #evergreen"
filename = f"Turbulence_{casename}-000000-000720.nc"
filename2 = f"EddyFlux_{casename}-000000-000720.nc"

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
    ds2 = xr.open_dataset(f"{DatasetDir}/{filename2}")
    wanothano = ds2['wanothano'][t0:t1+1, :, xIdx]
    
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


    return tke, enstrophy, wanothano, BLH1, BLH2, time_range, zc

tke0, enstrophy0, wanothano0, BLH10, BLH20, time_range, zc = loadData(t0, t1, xIdx=0)
tke1, enstrophy1, wanothano1, BLH11, BLH21, time_range, zc = loadData(t0, t1, xIdx=1)

#%% Draw BLH

#### USER ARGUMENT START #####

## NECESSARY ARGUMENT ##

# DATA
line_data = [BLH11, BLH21]
shade_data = wanothano1.T

## OPTIONAL ARGUMENT ##

# Fig Frame Setting
title = r"Vertical $\theta$ transport"
casename = casename
figsize = [6,4]
imgname = f"wth-{casename}.png"
savefig = True
show = False

# Axis Infromation Setting Time (x), Height (y)
zbnd = [0, 2000]
tbnd = [t0, t1]
toffset_min = 300
tresolution_min = 2

# Multiple or Single Line Data and Setting (May be single object or list)
line_color_list = ["limegreen","tab:orange"]
line_style_list = ["-", "-"]
line_width_list = [2.5, 1.8]
line_label_list = ["",""] #[r'$\theta_{sfc}+0.5K$', r'max $d\theta/dz$']

# Single Pcolormesh Data and Setting
shade_levels = np.arange(-0.1, 0.101, 0.01)
shade_cmap = 'coolwarm'
shade_extend = 'both'
shade_cbar_label = r"$\overline{w'\theta'}$ [K]"

#### USER ARGUMENT END #####

# Pre-processing
tstart=tbnd[0]
tend=tbnd[1]
tminute = toffset_min + tstart * tresolution_min
startHHMM = f"{str(tminute//60).zfill(2)}:{str(tminute%60).zfill(2)}"
time_range = pd.date_range(start=f"2024-01-01 {startHHMM}", periods=(tend-tstart+1), freq="2T")

# Plot
fig, ax  = plt.subplots(1,1, figsize=figsize, dpi=300)
ax.set_ylim(zbnd)
ax.set_ylabel('z [m]')
ax.set_xlabel('Local Time')
ax.set_xlim([time_range[0], time_range[-1]])
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.minorticks_on()

ax.set_title(f"{title}", loc="left", fontsize=14, weight="bold")
ax.set_title(f"{casename}", loc="right", fontsize=14, weight="bold")

# Plot Lines: BLH
for data, color, style, width, label in zip(line_data, line_color_list, 
                     line_style_list, line_width_list, line_label_list):
    ax.plot(time_range, data, color=color, ls=style, lw=width, label=label)
    
cmap = plt.cm.get_cmap(shade_cmap, len(shade_levels)+1)
norm = BoundaryNorm(shade_levels, ncolors=len(shade_levels)+1, clip=False, extend=shade_extend)
PC0=ax.pcolormesh(time_range, zc, shade_data, cmap=cmap, norm=norm)

# Plot Shaded: wth
[x0,y0],[x1,y1]=ax.get_position().get_points()
cax3 = fig.add_axes([x1+0.02, y0+0.002, 0.015, y1-y0-0.004]) # x0, y0, width, height
CB3 = plt.colorbar(PC0, cax=cax3,orientation='vertical', extend=shade_extend)
CB3.set_label(shade_cbar_label,labelpad=0.10, fontsize=12)     

# Add legend to plot
ax.legend(loc=2, ncol=1,fontsize=12,borderpad=0.2,handleheight=0.9,handlelength=1.5,handletextpad=0.3,
          labelspacing=0.2,columnspacing=1.0,framealpha=0.90)

if savefig: 
    plt.savefig(imgname, bbox_inches='tight')
if show:
    plt.show()
