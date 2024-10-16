#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 20:27:08 2024

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

#%%

DatasetDir = "/data/mlcloud/d11229002/VVM/DATA"
casename = "pbl_up"

#DatasetDir = "/data/yhc2080/VVM/DATA"
#casename = "pbl_half_PU_1" #evergreen_qc"
t0 = 0
t1 = 720
ntime = int(t1-t0+1)
nz = 50
ResultSlicerList = [slice(0,1),np.arange(50),np.arange(128),np.arange(128)]
TOPO = TaiwanVVMTOPO(f"{DatasetDir}/{casename}/TOPO.nc")



#%%

Landtype = TOPO.loadVariable('lu')

colors = [
    'purple',  # 1. Urban and built-up land
    '#FFD700',  # 2. Dryland, cropland, and pasture
    '#ADFF2F',  # 3. Irrigated cropland and pasture
    '#32CD32',  # 4. Mixed type of 2 and 3
    '#9ACD32',  # 5. Cropland/Grassland Mosaic
    '#556B2F',  # 6. Cropland/Woodland Mosaic
    '#7CFC00',  # 7. Grassland
    '#8B4513',  # 8. Shrubland
    '#A0522D',  # 9. Mixed type of 7 and 8
    '#D2691E',  # 10. Savanna
    '#228B22',  # 11. Deciduous broadleaf forest
    '#006400',  # 12. Deciduous needleleaf forest
    '#32CD32',  # 13. Evergreen broadleaf forest
    '#008000',  # 14. Evergreen needleleaf forest
    '#6B8E23',  # 15. Mixed forest
    '#4169E1',  # 16. Water bodies (ocean or lake)
    '#66CDAA',  # 17. Herbaceous wetland
    '#8FBC8F',  # 18. Wooded wetland
    '#D2B48C',  # 19. Barren and sparsely vegetated
    '#FFDAB9',  # 20. Herbaceous tundra
    '#CD853F',  # 21. Wooded tundra
    '#F5DEB3',  # 22. Mixed tundra
    '#FFE4B5',  # 23. Bare ground tundra
    '#FFFFFF'   # 24. Snow or ice
]

# 建立自定義的 colormap
cmap = ListedColormap(colors)

# 設定每一類別的邊界
boundaries = np.arange(1, 26)  # 因為有 24 種類型，邊界範圍 1 到 25
norm = BoundaryNorm(boundaries, cmap.N, clip=True)


fig, ax  = plt.subplots(1,1, figsize=[6,4], dpi=300, sharex=True, sharey=True)


# 使用 pcolormesh 繪製 Landtype 陣列
c = ax.pcolormesh(Landtype, cmap=cmap, norm=norm)

# 設定 colorbar
cbar = fig.colorbar(c, ax=ax, boundaries=boundaries, ticks=np.arange(1, 25)+0.5)
cbar.ax.set_yticklabels([f'{i}' for i in range(1, 25)])  # 設置 colorbar 標籤

# 顯示圖表
ax.set_title("Land Type", loc="left")
ax.set_title(casename, loc="right", weight ="bold")
plt.savefig(f"Landtype-{casename}.png")
plt.show()