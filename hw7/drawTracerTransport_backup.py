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
from matplotlib.ticker import MultipleLocator

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


xc = np.arange(0, 128*200, 200)+100

def generate_test_data(t0, t1, n_points=128):
    # 產生時間範圍 (t0 到 t1)，以 1 為間隔 (總共 t1-t0+1 個時間點)
    time_range = np.linspace(2*np.pi, 0, t1 - t0 + 1)

    # 產生空間範圍 (例如 128 個點)
    space_range = np.linspace(0, np.pi, n_points)

    # 使用 np.outer 來生成 2D array，基於 np.sin() 的時間和空間變化
    data = np.outer(np.sin(time_range), np.cos(space_range)) 
    
    return data

# 測試產生 (t1-t0+1, 128) 大小的 2D array
test_data = (generate_test_data(t0, t1, n_points=len(xc)) +1)*10

print(test_data.shape)  # 會顯示 (t1-t0+1, 128)


#%% Draw BLH

#### USER ARGUMENT START #####

## NECESSARY ARGUMENT ##

# DATA
shade_data = test_data
xc = np.arange(0, 128*200, 200)+100

## OPTIONAL ARGUMENT ##

# Fig Frame Setting
title = r"Vertical $\theta$ transport"
casename = casename
figsize = [6,4]
imgname = f"wth-{casename}.png"
savefig = True
show = False

# Axis Infromation Setting Time (x), Height (y)
xbnd = [0, 128*200]
tbnd = [t0, t1]
toffset_min = 300
tresolution_min = 2


# Single Pcolormesh Data and Setting
shade_levels = np.arange(2,20.1,2)
shade_cmap = 'YlOrBr'
shade_extend = 'both'
shade_cbar_label = r"Tracer [a.u.]"

#### USER ARGUMENT END #####

# Pre-processing
tstart=tbnd[0]
tend=tbnd[1]
tminute = toffset_min + tstart * tresolution_min
startHHMM = f"{str(tminute//60).zfill(2)}:{str(tminute%60).zfill(2)}"
time_range = pd.date_range(start=f"2024-01-01 {startHHMM}", periods=(tend-tstart+1), freq="2T")

# Plot
fig, ax  = plt.subplots(1,1, figsize=figsize, dpi=300)
ax.set_xlim(xbnd)
ax.xaxis.set_major_locator(MultipleLocator(3200))
ax.set_ylabel('Local Time')
ax.set_xlabel('x [m]')
ax.set_ylim([time_range[0], time_range[-1]])
ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.minorticks_on()

ax.set_title(f"{title}", loc="left", fontsize=14, weight="bold")
ax.set_title(f"{casename}", loc="right", fontsize=14, weight="bold")

    
cmap = plt.cm.get_cmap(shade_cmap, len(shade_levels)+1)
norm = BoundaryNorm(shade_levels, ncolors=len(shade_levels)+1, clip=False, extend=shade_extend)
PC0=ax.pcolormesh(xc, time_range, shade_data, cmap=cmap, norm=norm)

# Plot Shaded: Tracer
[x0,y0],[x1,y1]=ax.get_position().get_points()
cax3 = fig.add_axes([x1+0.02, y0+0.002, 0.015, y1-y0-0.004]) # x0, y0, width, height
CB3 = plt.colorbar(PC0, cax=cax3,orientation='vertical', extend=shade_extend)
CB3.set_label(shade_cbar_label,labelpad=0.10, fontsize=12)     

if savefig: 
    plt.savefig(imgname, bbox_inches='tight')
if show:
    plt.show()
