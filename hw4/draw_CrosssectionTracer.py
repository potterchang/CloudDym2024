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
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
#%% Load Functions
DatasetDir = "/data/yhc2080/VVM/DATA"
casename = "pbl_tracer"
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
    ind_y = 64

    u = VVMLoader.loadVariable("u", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    v = VVMLoader.loadVariable("v", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    w = VVMLoader.loadVariable("w", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)

    TKE = u**2+v**2+w**2    

    th = VVMLoader.loadVariable("th", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)

    rootgrp = nc.Dataset(f"{DatasetDir}/{casename}/archive/{casename}.L.Dynamic-{str(itime).zfill(6)}.nc")    
    eta = centerRegridder(rootgrp, 'eta', [slice(0,1),np.arange(50),ind_y,np.arange(128)])
    zeta = centerRegridder(rootgrp, 'zeta',[slice(0,1),np.arange(50),ind_y,np.arange(128)])
    xi = centerRegridder(rootgrp, 'xi', [slice(0,1),np.arange(50),ind_y,np.arange(128)])    
    Enstrophy=eta**2+zeta**2+xi**2
    
    tr01 = VVMLoader.loadVariable("tr01", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    tr02 = VVMLoader.loadVariable("tr02", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    tr03 = VVMLoader.loadVariable("tr03", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    
    return u, v, w, th, tr01, tr02, tr03, TKE, Enstrophy


u, v, w, th, tr01, tr02, tr03, TKE, Enstrophy = loadData(casename, itime)
#%%

def drawCross(casename, itime):
    toffset = 300
    tsteplength = 2
    tminute = toffset + itime * tsteplength
    HHMM = f"{str(tminute//60).zfill(2)}:{str(tminute%60).zfill(2)}"
    u, v, w, th, tr01, tr02, tr03, TKE, Enstrophy = loadData(casename, itime)
    xc = np.arange(0, 128*200, 200)
    yc = np.arange(0, 128*200, 200)
    zc = np.concatenate(([0], np.arange(20, 1940.1, 40)))
    trhigh = [500, 300, 20]
    trs = [tr03, tr02, tr01]
    
    fig, ax = plt.subplots(3,1, figsize=[8,6] , sharex=True,sharey=True, dpi=300)

    ax[0].set_title(f"{casename}\n y=12.8km Cross-section", loc="right")
    ax[0].set_title(f"{HHMM}", fontsize=16)
    for i in range(3):
        ax[i].set_title(f"Source Height: {trhigh[i]}m", fontsize=12, loc='left', weight='bold')
        ax[i].set_ylabel('zc [m]')
        
        levels = np.arange(2, 15.1, 1)#np.arange(0, 1.75, 0.2)
        cmap = plt.cm.get_cmap('YlOrBr', len(levels)+1)
        norm = BoundaryNorm(levels, ncolors=len(levels)+1, clip=False, extend='max')
        PC = ax[i].pcolormesh(xc, zc, trs[i], cmap=cmap, norm=norm)
        
        ax[i].contour(xc, zc, TKE, levels=[0.2], colors=['b'])
        #ax[i].contour(xc, zc, Enstrophy[0], levels=[2e-5], colors=['g']) 
        
        thbot = th[1, :] 
        P05K = thbot + 0.5
        ind_posi = np.array((th - P05K) >= 0)
        adj = np.sum(ind_posi, axis=0) == 0
        ind_BLH1 = np.argmax(ind_posi, axis=0) - adj
        ind_BLH1 = np.ma.array(ind_BLH1, mask=ind_BLH1<0)
        BLH1 = np.take(zc, ind_BLH1)
        ax[i].plot(xc, BLH1, 'g', lw=2.5)
        
    ax[i].set_xlim([8000, 17600])
    ax[-1].set_ylim([0,1500])
    ax[-1].set_xlabel('xc [m]')
    ax[-1].set_xticks(np.arange(8000,17600.1, 1600))
    [x0,y0],[x1,y1]=ax[-1].get_position().get_points()
    cax3 = fig.add_axes([x0, y0-0.09, x1-x0-0.008, 0.016]) # x0, y0, width, height
    CB3 = plt.colorbar(PC, cax=cax3,orientation='horizontal', extend='both')
    CB3.set_label(f'Tracer [unit]',labelpad=0.10)
    
    BLH1legend = Line2D([0], [0], linestyle='-', color='b',linewidth=1.5),
    BLH2legend = Line2D([0], [0], linestyle='-', color='g',linewidth=2.5),
    ax[0].legend([BLH1legend,BLH2legend],
                     [r'TKE 0.2 $m^2 s^{-2}$',r'$\theta_{sfc}$+0.5K'],loc=1,
                     fontsize=12,borderpad=0.3,handleheight=0.9,handlelength=1.5,handletextpad=0.4,
                     labelspacing=0.2,columnspacing=1.0,framealpha=0.90)
    
    foldername = f"crossTracer/{casename}"
    os.makedirs(foldername, exist_ok=True)
    imgname = f"crossTracer-{casename}-{str(itime).zfill(6)}.png"
    plt.savefig(f"{foldername}/{imgname}", bbox_inches="tight")
    
    return imgname

drawCross(casename, itime=240)


#%%
starttime = time.time()
if __name__ == '__main__':
    
    for casename in ["pbl_tracer"]:
        try:
            nProc = int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()/12)) # core to use
            with Pool(nProc) as p:
                results = [p.apply_async(drawCross, 
                    (casename, itime, )) for itime in range(t0,t1+1,1)]
                from tqdm import tqdm
                fin = [result.get() for result in tqdm(results)]
        except Exception as e:
            print("An error occurred:", e)
        print("finish", casename)
        print(f"Elapsed: {time.time()-starttime} sec.")

print("Done Everything.")