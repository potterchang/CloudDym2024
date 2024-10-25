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
nx = 128
ResultSlicerList = [slice(0,1),np.arange(50),np.arange(128),np.arange(128)]
TOPO = TaiwanVVMTOPO(f"{DatasetDir}/{casename}/TOPO.nc")
varlist = ["tr01","tr02","tr03","tr04","NO","NO2","O3"]
namelist= ["tracer01","tracer02","tracer03","tracer04","NO","NO2","O3"]
unitlist = ["a.u.","a.u.","a.u.","a.u.","ppb","ppb","ppb"]

def varMean(var_full):
    var = np.nanmean(var_full, axis=(1))
    return var

def calPollutants(casename, itime, varlist):
    
    VVMLoader = TaiwanVVMData(DatasetDir=DatasetDir, ExampleCasename=casename) 
    
    vardata = []
    for varname in varlist:
        var_full = VVMLoader.loadVariable(varname, itime, casename, zIdx=np.arange(50),yIdx=np.arange(128),xIdx=np.arange(128),TOPOmask=False, regrid=False)
        var = varMean(var_full)
        vardata.append(var)
    
    return vardata


def save_netcdf(vardata, save_path, casename, t0, t1):
    
    output_folder = f'{save_path}'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    filename = f"Pollutants_{casename}-{str(t0).zfill(6)}-{str(t1).zfill(6)}.nc"    
    
    with nc.Dataset(f"{output_folder}/{filename}", 'w') as ds:
        # Create dimensions
        if 'time' not in ds.dimensions:
            ds.createDimension('time', None)  # Unlimited dimension for time
        if 'z' not in ds.dimensions:
            ds.createDimension('z', 50)
        if 'x' not in ds.dimensions:
            ds.createDimension('x', 128)

        # Create and Set variables of coordinate
        time_var = ds.createVariable('time', 'i4', ('time'))
        time_var.setncatts({'name':'time','unit':'output timestep'})
        time_var[:] = np.arange(t0, t1+1)
        z_var = ds.createVariable('zc', 'i4', ('z'))
        z_var.setncatts({'name': "vertical height of model layers, MSL",'unit':'m'})
        z_var[:] = np.concatenate(([0], np.arange(20, 1940.1, 40)))
        x_var = ds.createVariable('xc', 'i4', ('x'))
        x_var.setncatts({'name': "x-coordinate of grid cell centers in Cartesian system",'unit':'m'})
        x_var[:] = np.arange(0, 128*200, 200) + 100
        
        for var, varname, varlongname, varunit in zip(vardata, varlist, namelist, unitlist):
            # Create variables        
            q_var = ds.createVariable(varname, 'f4', ('time', 'z', 'x'),
                                           compression='zlib', complevel=8)
            # Set attributes
            q_var.setncatts({'name':varlongname, 'unit':varunit})
            # Store data
            q_var[:] =var
        
        # Global attribute
        ds.setncatts({'comment': "y-axis averaged tracer and chemicals"})
       
    return filename

#%% Main

starttime = time.time()

if __name__ == '__main__':    

    try:
        nProc = 5 #int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()/12)) # core to use
        with Pool(nProc) as p:
            results=[p.apply_async(calPollutants,(casename, itime, varlist,)) for itime in range(t0, t1+1)]
            from tqdm import tqdm
            fin = [result.get() for result in tqdm(results)]
    except:
        print("Finish all")
        print(f"Elapsed: {time.time()-starttime} sec.")

VarsData = []
for i in range(len(varlist)):
    Vars = np.stack([vardata[i] for vardata in fin], axis=0)
    VarsData.append(Vars)

save_path = "../DATA"
filename = save_netcdf(VarsData, save_path, casename, t0, t1)
print("Done")
print(f"Elapsed: {time.time()-starttime} sec.")

