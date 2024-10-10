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
filename = f"Turbulence_{casename}-000000-000720.nc"

t0 = 60
t1 = 481
xIdx= 0
loc_list = ["Pasture","Urban"]
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

#%% Draw Variables

fig, ax  = plt.subplots(3,1, figsize=[8,9], dpi=300, sharex=True, sharey=True)

ax[0].set_ylim([0,2000])
ax[0].minorticks_on()
#ax[0].set_ylabel('z [m]')
ax[1].set_ylabel('z [m]')
ax[-1].set_xlabel('Local Time')
ax[-1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

ax[0].set_title(r'TKE [$m^{2} s^{-2}$]',loc='left')
ax[0].set_title(f"{casename}", loc="right", fontsize=14, weight="bold")
ax[0].set_title(f"{loc_list[xIdx]} Avg.", loc="center", fontsize=14, weight="bold")
levels = np.arange(0, 2.01, 0.2)#np.arange(0, 1.75, 0.2)
cmap = plt.cm.get_cmap('RdYlBu_r', len(levels)+1)
norm = BoundaryNorm(levels, ncolors=len(levels)+1, clip=False, extend='max')
PC0 = ax[0].pcolormesh(time_range, zc, tke.T, cmap=cmap, norm=norm)
plt.colorbar(PC0)

ax[1].set_title(r'Enstrophy [$10^{-5} s^{-2}$]',loc='left')
levels = np.arange(0, 11.01, 1.0)#np.arange(0, 8, 1)
cmap = plt.cm.get_cmap('RdYlBu_r', len(levels)+1)
norm = BoundaryNorm(levels, ncolors=len(levels)+1, clip=False, extend='max')
PC1 = ax[1].pcolormesh(time_range, zc, enstrophy.T*1e5, cmap=cmap, norm=norm)
plt.colorbar(PC1)

ax[2].set_title(r'$\theta$ [K]',loc='left')
levels = np.arange(294, 307, 1)#np.arange(290, 303, 1)
cmap = plt.cm.get_cmap('coolwarm', len(levels)+1)
norm = BoundaryNorm(levels, ncolors=len(levels)+1, clip=False, extend='both')
PC2 = ax[2].pcolormesh(time_range, zc, th.T, cmap=cmap, norm=norm)
plt.colorbar(PC2)

ax[2].plot(time_range, BLH1, color='limegreen', label=r'$\theta_{sfc}+0.5K$')
ax[2].plot(time_range, BLH2, color='k', label=r'max $d\theta/dz$')
ax[2].legend(loc=2)

plt.savefig(f"Evolution-{casename}-{loc_list[xIdx]}.png", bbox_inches='tight')

