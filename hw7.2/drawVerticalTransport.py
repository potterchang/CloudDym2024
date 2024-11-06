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
sys.path.append("/data/yhc2080/UTIL/VVMAnalysisKit")
from plotter import PBLplot
# Visualization
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
#%% Load Data
DatasetDir = "../DATA"
casename = "pbl_half_PU_uarea_2" #evergreen"
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

    ## BLH3: eddy 0.01 K
    BLH3 = np.full(th.shape[0], np.nan)
    # 從底層往上查找溫度梯度變化
    for i in range(wanothano.shape[0]):  # 對於每一個時間步
        if wanothano[i, 0] >= 0.001:  # 底層溫度 >= 0.01K
            # 找到從底層往上第一次小於 0.01K 的高度層
            index = np.where(wanothano[i, :] < 0.001)[0]
            if index.size > 0:  # 如果找到了小於 0.01K 的層
                BLH3[i] = zc[index[0]]  # 對應的高度
        else:
            BLH3[i] = np.nan  # 如果底層溫度 < 0.01K，設為 np.nan


    return tke, enstrophy, wanothano, BLH1, BLH2, BLH3, time_range, zc

tke0, enstrophy0, wanothano0, BLH10, BLH20, BLH30, time_range, zc = loadData(t0, t1, xIdx=0)
tke1, enstrophy1, wanothano1, BLH11, BLH21, BLH31, time_range, zc = loadData(t0, t1, xIdx=1)

#%% Draw BLH
xIdx = 1
PBLheight = [BLH11,BLH21]
wth = wanothano1.T
imgname = f'Wth-{casename}-{loc_list[xIdx]}.png'

PBLplot().plot_vertical_theta_transport(PBLheight, wth, zc, 
                                  casename=f"{casename[9:]}", title=f"{loc_list[xIdx]}", 
                                  figsize=[5,4], imgname=imgname, savefig=True, show=True, 
                                  zbnd=[0, 2000], 
                                  tbnd=[t0, t1], toffset_min=300, tresolution_min=2, 
                                  line_color_list=None, line_style_list="-", 
                                  line_width_list=2.5, line_label_list=[r'$\theta_{sfc}+0.5K$', r'max $d\theta/dz$'], 
                                  shade_levels=np.arange(-0.1, 0.101, 0.01))
