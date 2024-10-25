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
#casename="pbl_up"

DatasetDir = "/data/yhc2080/VVM/DATA"
casename = "pbl_half_PU_uarea_2" #evergreen_qc"
t0 = 0
t1 = 720
ntime = int(t1-t0+1)
nz = 50
ResultSlicerList = [slice(0,1),np.arange(50),np.arange(128),np.arange(128)]
TOPO = TaiwanVVMTOPO(f"{DatasetDir}/{casename}/TOPO.nc")

def varMean(var_full):
    var1 = np.nanmean(var_full[:,0:64], axis=(0,1))
    var2 = np.nanmean(var_full[:,64:128], axis=(0,1))
    var = np.stack([var1, var2], axis=0)
    return var

def calEnergy(casename, itime):
    
    VVMLoader = TaiwanVVMData(DatasetDir=DatasetDir, ExampleCasename=casename)    
    wth_full = VVMLoader.loadVariable("wth", itime, casename, zIdx=np.arange(50),yIdx=np.arange(128),xIdx=np.arange(128),TOPOmask=False, regrid=False)
    wqv_full = VVMLoader.loadVariable("wqv", itime, casename, zIdx=np.arange(50),yIdx=np.arange(128),xIdx=np.arange(128),TOPOmask=False, regrid=False)
    fdswtoa_full = VVMLoader.loadVariable("fdswtoa", itime, casename, zIdx=np.arange(50),yIdx=np.arange(128),xIdx=np.arange(128),TOPOmask=False, regrid=False)
 
    wth = varMean(wth_full)
    wqv = varMean(wqv_full)
    fdswtoa = varMean(fdswtoa_full)
    
    return wth, wqv, fdswtoa


def save_netcdf(wths, wqvs, fdswtoas, save_path, casename, t0, t1):
    
    output_folder = f'{save_path}'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    filename = f"Energy_{casename}-{str(t0).zfill(6)}-{str(t1).zfill(6)}.nc"    
    
    with nc.Dataset(f"{output_folder}/{filename}", 'w') as ds:
        # Create dimensions
        if 'time' not in ds.dimensions:
            ds.createDimension('time', None)  # Unlimited dimension for time
        if 'x' not in ds.dimensions:
            ds.createDimension('x', 2)

        # Create variables
        time_var = ds.createVariable('time', 'i4', ('time'))
        wth_var = ds.createVariable('wth', 'f4', ('time', 'x'),
                                       compression='zlib', complevel=8)
        wqv_var = ds.createVariable('wqv', 'f4', ('time', 'x'),
                                       compression='zlib', complevel=8)
        fdswtoa_var = ds.createVariable('fdswtoa', 'f4', ('time', 'x'),
                                       compression='zlib', complevel=8)

        # Set attributes
        time_var.setncatts({'name':'time','unit':'output timestep'})
        wth_var.setncatts({'name':'Surface flux of potential temperature', 'unit':'K kg m-2 s-1'})
        wqv_var.setncatts({'name':'Surface flux of water vapor', 'unit':'kg m-2 s-1'})        
        fdswtoa_var.setncatts({'name':'Downward flux of shortwave radiation at TOA', 'unit':'W/m**2'})  

        # Store data
        time_var[:] = np.arange(t0, t1+1)
        wth_var[:] = wths
        wqv_var[:] = wqvs
        fdswtoa_var[:] = fdswtoas

        
        # Global attribute
        ds.setncatts({'comment': "regional averaged energy budget"})
       
    return filename

#%% Main

starttime = time.time()

if __name__ == '__main__':    

    try:
        nProc = 5 #int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()/12)) # core to use
        with Pool(nProc) as p:
            results=[p.apply_async(calEnergy,(casename, itime,)) for itime in range(t0, t1+1)]
            from tqdm import tqdm
            fin = [result.get() for result in tqdm(results)]
    except:
        print("Finish all")
        print(f"Elapsed: {time.time()-starttime} sec.")

wths = np.full([ntime, 2], np.nan) 
wqvs = np.full([ntime, 2], np.nan) 
fdswtoas = np.full([ntime, 2], np.nan) 

for itime in range(ntime):

    wths[itime] = fin[itime][0]
    wqvs[itime] = fin[itime][1]
    fdswtoas[itime] = fin[itime][2]

save_path = "../DATA"
filename = save_netcdf(wths, wqvs, fdswtoas, save_path, casename, t0, t1)
print("Done")
print(f"Elapsed: {time.time()-starttime} sec.")

