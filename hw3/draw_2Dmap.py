#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 18:40:19 2024

@author: yhc2080
"""

#%% Initialize
# Basic
import os
import sys
import time
import numpy as np
#import geopandas as gpd
# I/O Processing
import pickle
from scipy.interpolate import interpn, RegularGridInterpolator
import netCDF4 as nc
import multiprocessing
from multiprocessing import Pool
# Taiwan VVM
sys.path.append("/data/yhc2080/UTIL")
from TaiwanVVMLoader import TaiwanVVMTOPO, TaiwanVVMData
# Visualization
import matplotlib.pyplot as plt

#%% Load Functions
DatasetDir = "/data/yhc2080/VVM/DATA"
casename = "pbl_evergreen_qc"
t0 = 0
t1 = 420
ntime = int(t1-t0+1)
nz = 50
ResultSlicerList = [slice(0,1),np.arange(50),np.arange(128),np.arange(128)]
TOPO = TaiwanVVMTOPO(f"{DatasetDir}/{casename}/TOPO.nc")
    
def centerRegridder(rootgrp, varName, ResultSlicerList, 
                    ResultDimList=None, zbot=False, #for new zbot method only
                    lat_bnd=None, lon_bnd=None, yIdx=None, xIdx=None, #for new zbot method only
                    TOPOmask=True, TOPOmethod='nearest', fill_value=None, keepdims=False):

    # If varName is not in Famliy of displaced variables, 
    # return the identical array.
    
    ### In (z,y,x) dimension ###
    # Family I: Velocity field
    mask_roll = {'u':-1, 'v':-1, 'w':-1}
    cal_roll  = {'u': 1, 'v': 1, 'w': 1}
    axis_roll = {'u': 2, 'v': 1, 'w': 0}
    
    # Family II: Vorticity field
    mask_roll1 = {'xi':-1, 'eta':-1, 'zeta':-1}
    cal_roll1  = {'xi': 1, 'eta': 1, 'zeta': 1}
    axis_roll1 = {'xi': 0, 'eta': 0, 'zeta': 1}
    mask_roll2 = {'xi':-1, 'eta':-1, 'zeta':-1}
    cal_roll2  = {'xi': 1, 'eta': 1, 'zeta': 1}
    axis_roll2 = {'xi': 1, 'eta': 2, 'zeta': 2}
    
    # TOPO mask value setting
    mask_values = {'nearest': np.nan, 'linear0': 0}
    if fill_value is None:
        fill_value = mask_values[TOPOmethod]
    
    # Regridding
    Var_ori_sliced = rootgrp.variables[varName][tuple(ResultSlicerList)]
    varshape = rootgrp.variables[varName].shape
    
    if varName in ['u','v','w']:
        ### Single Rolled ###     
        # Create rolled slicer
        idim = axis_roll[varName] +1
        if type(ResultSlicerList[idim])==slice:
            dim_indices = list(range(ResultSlicerList[idim].start, ResultSlicerList[idim].stop))
            dim_indices = [(idx - cal_roll[varName]) % varshape[idim] for idx in dim_indices]
        else:
            dim_indices = (ResultSlicerList[idim] - cal_roll[varName]) % varshape[idim]
        RolledSlicerList = list(ResultSlicerList)
        RolledSlicerList[idim] = dim_indices
        RolledSlicerList = tuple(RolledSlicerList)
        # Load rolled data
        Var_rolled_sliced = rootgrp.variables[varName][tuple(RolledSlicerList)]

        # Calculate regridded value
        OUT = np.nanmean([Var_ori_sliced,Var_rolled_sliced],axis=0)    
    elif varName in ['xi','eta','zeta']:
        ### Triple Rolled ###
        # Create rolled slicer
        idim1 = axis_roll1[varName] +1
        if type(ResultSlicerList[idim1])==slice:
            dim1_indices = list(range(ResultSlicerList[idim1].start, ResultSlicerList[idim1].stop))
            dim1_indices = [(idx - cal_roll1[varName]) % varshape[idim1] for idx in dim1_indices]
        else:
            dim1_indices = (ResultSlicerList[idim1] - cal_roll1[varName]) % varshape[idim1]
        idim2 = axis_roll2[varName] +1
        if type(ResultSlicerList[idim2])==slice:
            dim2_indices = list(range(ResultSlicerList[idim2].start, ResultSlicerList[idim2].stop))
            dim2_indices = [(idx - cal_roll2[varName]) % varshape[idim2] for idx in dim2_indices]
        else:
            dim2_indices = (ResultSlicerList[idim2] - cal_roll2[varName]) % varshape[idim2]
        RolledSlicerList1 = list(ResultSlicerList)
        RolledSlicerList1[idim1] = dim1_indices
        RolledSlicerList1 = tuple(RolledSlicerList1)
        RolledSlicerList2 = list(ResultSlicerList)
        RolledSlicerList2[idim2] = dim2_indices
        RolledSlicerList2 = tuple(RolledSlicerList2)
        RolledSlicerList12 = list(ResultSlicerList)
        RolledSlicerList12[idim1] = dim1_indices
        RolledSlicerList12[idim2] = dim2_indices
        RolledSlicerList12 = tuple(RolledSlicerList12)
        # Load rolled data
        Var_rolled_sliced1 = rootgrp.variables[varName][tuple(RolledSlicerList1)]
        Var_rolled_sliced2 = rootgrp.variables[varName][tuple(RolledSlicerList2)]
        Var_rolled_sliced12 = rootgrp.variables[varName][tuple(RolledSlicerList12)]

        # Calculate regridded value
        OUT = np.nanmean([Var_ori_sliced,Var_rolled_sliced1,Var_rolled_sliced2,Var_rolled_sliced12],axis=0)    
    
    else:
        OUT = Var_ori_sliced

    return OUT

def loadData(casename, itime):
    
    VVMLoader = TaiwanVVMData(DatasetDir=DatasetDir, ExampleCasename=casename)    
    zc = np.concatenate(([0], np.arange(20, 1940.1, 40)))
    ind_z = np.argmin(abs(zc-1000))    

    u = VVMLoader.loadVariable("u", itime, casename, zIdx=ind_z, yIdx=np.arange(128), xIdx=np.arange(128), TOPOmask=False, regrid=True)
    v = VVMLoader.loadVariable("v", itime, casename, zIdx=ind_z, yIdx=np.arange(128), xIdx=np.arange(128), TOPOmask=False, regrid=True)
    w = VVMLoader.loadVariable("w", itime, casename, zIdx=ind_z, yIdx=np.arange(128), xIdx=np.arange(128), TOPOmask=False, regrid=True)
 
    qv = VVMLoader.loadVariable("qv", itime, casename, zIdx=ind_z, yIdx=np.arange(128), xIdx=np.arange(128), TOPOmask=False, regrid=True)
    th = VVMLoader.loadVariable("th", itime, casename, zIdx=ind_z, yIdx=np.arange(128), xIdx=np.arange(128), TOPOmask=False, regrid=True)
    PIBAR = VVMLoader.loadFort98('PIBAR')
    t = th * PIBAR[ind_z]
    return u, v, w, qv, t

#%%

def drawMap2d(casename, itime):
    toffset = 300
    tsteplength = 2
    tminute = toffset + itime * tsteplength
    HHMM = f"{str(tminute//60).zfill(2)}:{str(tminute%60).zfill(2)}"
    u, v, w, qv, t = loadData(casename, itime)
    xc = np.arange(0, 128*200, 200)
    yc = np.arange(0, 128*200, 200)
    
    
    fig, ax = plt.subplots(1,2, figsize=[10,5] , sharex=True,sharey=True, dpi=300)
    ax[0].set_aspect(1)
    ax[1].set_aspect(1)
    
    CF1=ax[0].contourf(xc, yc, t, cmap='turbo', levels=np.arange(293.0, 294.51, 0.1), extend='both')
    plt.colorbar(CF1, label='T [K]', orientation='horizontal', extend='both')
    CT=ax[0].contour(xc, yc, qv*1000, levels=np.arange(0, 21.1, 3), colors="w", linewidths=0.8)
    plt.clabel(CT, fontsize=8)
    
    xygap = 5
    QV=ax[1].quiver(xc[::xygap], yc[::xygap], u[::xygap,::xygap], v[::xygap,::xygap], color="k", scale=15,zorder=10)
    ax[1].quiverkey(QV, 1.02, 1.00, 1.0, "\n1.0"+"\nm/s",angle=90, color="k",
                 labelpos = "E", labelcolor="k", fontproperties={'size':8.5}, zorder = 20)
    CF2=ax[1].contourf(xc, yc, w, cmap='turbo', levels=np.arange(-0.5, 0.51, 0.1), extend='both',zorder=5)
    plt.colorbar(CF2, label='w [m/s]', orientation='horizontal', extend='both')
    
    ax[0].set_title(f"                                        {casename}      {HHMM}"+"\nT (shade); QV (contour, 3g/kg interval)", loc="left")
    ax[1].set_title(f"\nW (shade); U, V (vector)", loc="left")
    ax[0].set_xticks(np.arange(0,25600.1, 6400))
    ax[0].set_yticks(np.arange(0,25600.1, 6400))
    
    foldername = f"map2d/{casename}"
    os.makedirs(foldername, exist_ok=True)
    imgname = f"map2d-{casename}-{str(itime).zfill(6)}.png"
    plt.savefig(f"{foldername}/{imgname}", bbox_inches="tight")
    return imgname

#%%
starttime = time.time()
if __name__ == '__main__':
    
    for casename in ["pbl_ctl", "pbl_evergreen_qc", "pbl_evergreen"]:
        try:
            nProc = int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()/12)) # core to use
            with Pool(nProc) as p:
                results = [p.apply_async(drawMap2d, 
                    (casename, itime, )) for itime in range(t0,t1+1,1)]
                from tqdm import tqdm
                fin = [result.get() for result in tqdm(results)]
        except Exception as e:
            print("An error occurred:", e)
        print("finish", casename)
        print(f"Elapsed: {time.time()-starttime} sec.")

print("Done Everything.")