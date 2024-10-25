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
from matplotlib.ticker import MultipleLocator
#%% Load Functions
DatasetDir = "/data/yhc2080/VVM/DATA"
casename = "pbl_half_PU_uarea_1"
t0 = 0
t1 = 720
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
    ind_y = np.arange(0,128)

    u = VVMLoader.loadVariable("u", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    v = VVMLoader.loadVariable("v", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    w = VVMLoader.loadVariable("w", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)

    TKE = u**2+v**2+w**2  

    u = np.nanmean(u, axis=1)
    v = np.nanmean(v, axis=1)
    w = np.nanmean(w, axis=1)
    TKE = np.nanmean(TKE, axis=1)
      

    th = VVMLoader.loadVariable("th", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    th = np.nanmean(th, axis=1)

    rootgrp = nc.Dataset(f"{DatasetDir}/{casename}/archive/{casename}.L.Dynamic-{str(itime).zfill(6)}.nc")    
    eta = centerRegridder(rootgrp, 'eta', [slice(0,1),np.arange(50),ind_y,np.arange(128)])
    zeta = centerRegridder(rootgrp, 'zeta',[slice(0,1),np.arange(50),ind_y,np.arange(128)])
    xi = centerRegridder(rootgrp, 'xi', [slice(0,1),np.arange(50),ind_y,np.arange(128)])    
    Enstrophy=eta**2+zeta**2+xi**2
    Enstrophy=np.nanmean(Enstrophy, axis=2)
    
    tr01 = VVMLoader.loadVariable("tr01", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    O3 = VVMLoader.loadVariable("O3", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    NO2 = VVMLoader.loadVariable("NO2", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)

    tr01 = np.nanmean(tr01, axis=1)
    O3 = np.nanmean(O3, axis=1)
    NO2 = np.nanmean(NO2, axis=1)

    #tr02 = VVMLoader.loadVariable("tr02", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    #tr03 = VVMLoader.loadVariable("tr03", itime, casename, zIdx=np.arange(50), yIdx=ind_y, xIdx=np.arange(128), TOPOmask=False, regrid=True)
    
    return u, v, w, th, tr01, O3, NO2, TKE, Enstrophy


u, v, w, th, tr01, O3, NO2, TKE, Enstrophy = loadData(casename, itime=90)
#%%

def drawCross(casename, itime):
    toffset = 300
    tsteplength = 2
    tminute = toffset + itime * tsteplength
    HHMM = f"{str(tminute//60).zfill(2)}:{str(tminute%60).zfill(2)}"
    u, v, w, th, tr01, O3, NO2, TKE, Enstrophy = loadData(casename, itime)
    xc = np.arange(0, 128*200, 200)+100
    yc = np.arange(0, 128*200, 200)
    zc = np.concatenate(([0], np.arange(20, 1940.1, 40)))
    trs = [tr01, NO2, O3]
    names = ["Tracer","NO2","O3"]
    levels_list = [np.arange(1.0,5.51,0.5)*200, np.arange(2,20.1,2)*100, np.arange(5,50.1,5)]
    cmap_list = ['YlOrBr','Purples','Blues']
    unit_list = ['a.u.','ppb','ppb']


    thbot = th[1, :] 
    P05K = thbot + 0.5
    ind_posi = np.array((th - P05K) >= 0)
    adj = np.sum(ind_posi, axis=0) == 0
    ind_BLH1 = np.argmax(ind_posi, axis=0) - adj
    ind_BLH1 = np.ma.array(ind_BLH1, mask=ind_BLH1<0)
    BLH1 = np.take(zc, ind_BLH1)
    
    fig, ax = plt.subplots(3,1, figsize=[7,7] , sharex=True,sharey=True, dpi=300)
    plt.subplots_adjust(hspace=0.25)
    ax[0].set_title(f"{casename}\n y-axis avg. Cross-section", loc="right",fontsize=10)
    ax[0].set_title(f"{HHMM}", fontsize=16)
    for i in range(3):
        ax[i].set_title(f"{names[i]} [{unit_list[i]}]", fontsize=14, loc='left', weight='bold')
        ax[i].set_ylabel('zc [m]')
        
        levels = levels_list[i]
        cmap = plt.cm.get_cmap(cmap_list[i], len(levels)+1)
        norm = BoundaryNorm(levels, ncolors=len(levels)+1, clip=False, extend='max')
        PC = ax[i].pcolormesh(xc, zc, trs[i], cmap=cmap, norm=norm, shading='nearest')
        
        xgap = 5; zgap=4
        wsmask = (u**2+w**2)**0.5 <= 0.1
        um = u.copy(); wm = w.copy()
        um[wsmask] = np.nan
        wm[wsmask] = np.nan
        QV=ax[i].quiver(xc[::xgap], zc[::zgap], um[::zgap,::xgap], wm[::zgap,::xgap], scale=40)
        if i ==0:
            QVK=ax[i].quiverkey(QV, 1.05, 1.1, 1, '1 m/s', coordinates='axes', 
                                labelsep=0.08, fontproperties={'size':9.5})

        ax[i].plot(xc, BLH1, 'grey', lw=2.5)
        ax[i].axvline(x=12800,color='r')

        [x0,y0],[x1,y1]=ax[i].get_position().get_points()
        cax3 = fig.add_axes([x1+0.02, y0+0.002, 0.012, y1-y0-0.004]) # x0, y0, width, height
        CB3 = plt.colorbar(PC, cax=cax3,orientation='vertical', extend='both')
        #CB3.set_label(f'[{unit_list[i]}]',labelpad=0.10)     
        ax[i].minorticks_on()
        
    ax[i].set_xlim([0, 25600])
    ax[-1].xaxis.set_major_locator(MultipleLocator(3200))
    ax[-1].set_ylim([0,1800])
    ax[-1].set_xlabel('xc [m]')

    BLH2legend = Line2D([0], [0], linestyle='-', color='grey',linewidth=2.5),
    ax[0].legend([BLH2legend],
                     [r'$\theta_{sfc}$+0.5K'],loc=2,
                     fontsize=12,borderpad=0.3,handleheight=0.9,handlelength=1.5,handletextpad=0.4,
                     labelspacing=0.2,columnspacing=1.0,framealpha=0.90)
    
    foldername = f"crossTracer/{casename}/Y64pm16"
    os.makedirs(foldername, exist_ok=True)
    imgname = f"crossTracer-{casename}-{str(itime).zfill(6)}.png"
    #plt.savefig(f"{foldername}/{imgname}", bbox_inches="tight")
    plt.show()
    return imgname

#drawCross(casename, itime=220)
drawCross("pbl_half_PU_uarea_1", itime=560)

#%%
starttime = time.time()
if __name__ == '__main__':
    
    for casename in ["pbl_half_PU_1","pbl_half_PU_2"]:
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
