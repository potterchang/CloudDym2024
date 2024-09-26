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
casename = "pbl_evergreen"
filename = f"Turbulence_{casename}-000000-000420.nc"

ds = xr.open_dataset(f"{DatasetDir}/{filename}")
time = ds['time']
zc = ds['zc']
tke = ds['tke'][0:421]
enstrophy = ds['enstrophy'][0:421]
th = ds['th'][0:421]
time_range = pd.date_range(start="2024-01-01 05:00", end="2024-01-01 19:00", freq="2T")

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
levels = np.arange(0, 1.1, 0.1)#np.arange(0, 1.75, 0.2)
cmap = plt.cm.get_cmap('RdYlBu_r', len(levels)+1)
norm = BoundaryNorm(levels, ncolors=len(levels)+1, clip=False, extend='max')
PC0 = ax[0].pcolormesh(time_range, zc, tke.T, cmap=cmap, norm=norm)
plt.colorbar(PC0)

ax[1].set_title(r'Enstrophy [$10^{-5} s^{-2}$]',loc='left')
levels = np.arange(0, 5.51, 0.5)#np.arange(0, 8, 1)
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
ax[2].legend()

plt.savefig(f"Evolution2-{casename}.png", bbox_inches='tight')

#%% Draw BLH

fig, ax  = plt.subplots(1,1, figsize=[6,4], dpi=300, sharex=True, sharey=True)

ax.set_ylim([0,1300])
ax.set_ylabel('z [m]')
ax.set_xlabel('Local Time')
ax.set_xlim([time_range[90], time_range[400]])
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

ax.set_title(f"{casename}", loc="right", fontsize=14, weight="bold")
ax.set_title(r'PBL Depth',loc='left')
levels = np.arange(0, 1.1, 0.2)
cmap = plt.cm.get_cmap('hot', len(levels)+1)
norm = BoundaryNorm(levels, ncolors=len(levels)+1, clip=False, extend='max')
PC0 = ax.contour(time_range, zc, tke.T, cmap=cmap, norm=norm)
ax.clabel(PC0, PC0.levels, inline=True, fmt='%1.1f', fontsize=12)
PC0_legend = Line2D([0], [0], linestyle='-', color='tab:red',linewidth=2.0)

levels = np.arange(0, 7, 1)
cmap = plt.cm.get_cmap('Blues', len(levels)+1)
norm = BoundaryNorm(levels, ncolors=len(levels)+1, clip=False, extend='max')
PC1 = ax.contour(time_range, zc, enstrophy.T*1e5, cmap=cmap, norm=norm)
ax.clabel(PC1, PC1.levels, inline=True, fmt='%1.0f', fontsize=12)
PC1_legend = Line2D([0], [0], linestyle='-', color='tab:blue',linewidth=2.0)


L1=ax.plot(time_range, BLH1, color='limegreen', lw=2.5, label=r'$\theta_{sfc}+0.5K$')
L2=ax.plot(time_range, BLH2, color='magenta', lw=2.5, label=r'max $d\theta/dz$')

# Combine all legend elements
all_handles = [PC0_legend, PC1_legend, L1[0], L2[0]]
all_labels = ['TKE', 'Enstrophy', r'$\theta_{sfc}+0.5K$', 'max $d\theta/dz$']

# Add legend to plot
ax.legend(all_handles, all_labels)

#ax.legend([PC0_legend, PC1_legend], ['TKE','Enstrophy'])

plt.savefig(f"PBL2-{casename}.png", bbox_inches='tight')
