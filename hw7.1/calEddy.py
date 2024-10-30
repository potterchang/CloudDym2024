#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:37:55 2024

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

#%% Load Functions
#DatasetDir = "/data/mlcloud/d11229002/VVM/DATA"
#casename = "pbl_up"

DatasetDir = "/data/yhc2080/VVM/DATA"
casename = "pbl_half_PU_uarea_1" #evergreen_qc"
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
 

def calEddyFlux(casename, itime):
    
    VVMLoader = TaiwanVVMData(DatasetDir=DatasetDir, ExampleCasename=casename)    
    ResultSlicerList = [slice(0,1),np.arange(50),np.arange(128),np.arange(128)]
    
    #u = VVMLoader.loadVariable("u", itime, casename, zIdx=np.arange(50), yIdx=np.arange(128), xIdx=np.arange(128), TOPOmask=False, regrid=True)
    #v = VVMLoader.loadVariable("v", itime, casename, zIdx=np.arange(50), yIdx=np.arange(128), xIdx=np.arange(128), TOPOmask=False, regrid=True)
    w = VVMLoader.loadVariable("w", itime, casename, zIdx=np.arange(50), yIdx=np.arange(128), xIdx=np.arange(128), TOPOmask=False, regrid=True)
    th = VVMLoader.loadVariable("th", itime, casename, zIdx=np.arange(50), yIdx=np.arange(128), xIdx=np.arange(128), TOPOmask=False, regrid=True)

    wbar = np.nanmean(w, axis=(1,2))
    wano = w - wbar.reshape(-1,1,1)
    
    thbar = np.nanmean(th, axis=(1,2))
    thano = th - thbar.reshape(-1,1,1)    
    
    EddyFlux = wano * thano
    
    Eddy1 = np.nanmean(EddyFlux[:,:,0:64], axis=(1,2))
    Eddy2 = np.nanmean(EddyFlux[:,:,64:128], axis=(1,2))
    Eddyfull = np.nanmean(EddyFlux[:,:,:], axis=(1,2))
    Eddy = np.stack([Eddy1, Eddy2, Eddyfull], axis=1)
    
    return Eddy


def save_netcdf(Eddys, save_path, casename, t0, t1):
    
    output_folder = f'{save_path}'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    filename = f"EddyFlux_{casename}-{str(t0).zfill(6)}-{str(t1).zfill(6)}.nc"    
    
    with nc.Dataset(f"{output_folder}/{filename}", 'w') as ds:
        # Create dimensions
        if 'time' not in ds.dimensions:
            ds.createDimension('time', None)  # Unlimited dimension for time
        if 'z' not in ds.dimensions:
            ds.createDimension('z', Eddys.shape[1])
        if 'x' not in ds.dimensions:
            ds.createDimension('x', Eddys.shape[2])

        # Create variables
        time_var = ds.createVariable('time', 'i4', ('time'))
        z_var = ds.createVariable('zc', 'i4', ('z'))
        Eddy_var = ds.createVariable('wanothano', 'f4', ('time', 'z', 'x'),
                                       compression='zlib', complevel=8)

        # Set attributes
        time_var.setncatts({'name':'time','unit':'output timestep'})
        z_var.setncatts({'name': "vertical height of model layers, MSL",'unit':'m'})
        Eddy_var.setncatts({'name':'eddy flux - w anomaly x th anomaly', 'unit':'K m s-1'})

        # Store data
        time_var[:] = np.arange(721)
        z_var[:] = np.concatenate(([0], np.arange(20, 1940.1, 40)))
        Eddy_var[:] = Eddys
        
        # Global attribute
        ds.setncatts({'comment': "3 averaged vertical structure: region 1 (x=1-64), region 2 (x=65-128), full domain"})
       
    return filename

#%% Main

starttime = time.time()

if __name__ == '__main__':    

    try:
        nProc = 5 #int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()/12)) # core to use
        with Pool(nProc) as p:
            results=[p.apply_async(calEddyFlux,(casename, itime,)) for itime in range(ntime)]
            from tqdm import tqdm
            fin = [result.get() for result in tqdm(results)]
    except:
        print("Finish all")
        print(f"Elapsed: {time.time()-starttime} sec.")

Eddys = np.full([ntime, nz, 3], np.nan) 

for itime in range(ntime):

    Eddys[itime] = fin[itime]

save_path = "../DATA"
filename = save_netcdf(Eddys, save_path, casename, t0, t1)
print("Done")
print(f"Elapsed: {time.time()-starttime} sec.")
