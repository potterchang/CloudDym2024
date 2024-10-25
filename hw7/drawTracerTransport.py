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
from matplotlib.ticker import MultipleLocator

#%% Load Data
DatasetDir = "../DATA"
casename = "pbl_half_PU_uarea_1" #evergreen"

t0 = 30
t1 = 720
z_interest = 500
poll_name_list = ['tr01', 'NO','NO2','O3']
cbar_label_list = ['Tracer [a.u.]','NO [ppb]','NO2 [ppb]','O3 [ppb]']
cmap_list = ['YlOrBr','Purples','Purples','Blues']
levels_list = [np.arange(100,1500.1,100), np.arange(10,100.1,5)*5, np.arange(10,100.1,5)*5, np.arange(5,60.1,5)*5]


def main(casename, poll_name, cbar_label, cmap_name, levels):

    filename = f"Pollutants_{casename}-000000-000720.nc"
    ds = xr.open_dataset(f"{DatasetDir}/{filename}")
    time = ds['time']
    zc = ds['zc']
    xc = ds['xc']
    imgname = f'Hov-{casename}-{poll_name}-zc{z_interest}.png'
    
    indz = np.argmin(np.array(abs(zc-z_interest)))
    pollutant = ds[poll_name][t0:t1+1, indz, :]
    PBLplot().plot_hovmoller_pollutants(pollutant, xc, 
                              casename=casename[9:], title=f"z= {z_interest} m", 
                              figsize=[3.6,5.4], imgname=imgname, savefig=True, show=True,
                              xbnd=[0, 25600], xLocator=6400,
                              tbnd=[t0,t1], toffset_min=300, tresolution_min=2, 
                              shade_levels=levels, 
                              shade_cmap=cmap_name, shade_extend='both', 
                              shade_cbar_label=cbar_label)

for poll_name, cbar_label, cmap_name, levels in zip(poll_name_list,cbar_label_list,cmap_list, levels_list):
    main(casename, poll_name, cbar_label, cmap_name, levels)   
