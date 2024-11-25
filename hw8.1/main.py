import numpy as np
import xarray as xr
#from plottools import dataPlotters
import matplotlib.pyplot as plt

# TODO 1: import your vvmtools
#from vvmtools import VVMTools
import vvmtools
# TODO 2: change the expname to your experiment name
# prepare expname and data coordinate
MotherDir = "/data/yhc2080/VVM/DATA"
expname  = 'pbl_PU_s1'
#expname = "pbl_half_PU_uarea_2"
#MotherDir = "/data/chung0823/VVM_cloud_dynamics_2024/DATA"
#expname = "pbl_hetero_dthdz_8"
nx = 128; x = np.arange(nx)*0.2
ny = 128; y = np.arange(ny)*0.2
nz = 50;  z = np.arange(nz)*0.04
nt = 721; t = np.arange(nt)*np.timedelta64(2,'m')+np.datetime64('2024-01-01 05:00:00')

# TODO 3: change the data to your data (seven lines for BL height and one shading for w'th')

# Read regional averaged data
DatasetDir = "../DATA"
filename = f"Turbulence_{expname}-000000-000720.nc"
filename2 = f"EddyFlux_{expname}-000000-000720.nc"

t0 = 0
t1 = 721
xIdx= 2
loc_list = ["pasture","urban", "domain"]

#vtls = vvmtools.VVMTools(case_path=f"{MotherDir}/{expname}")

def loadData(t0, t1, xIdx):
    ds = xr.open_dataset(f"{DatasetDir}/{filename}")
    time = ds['time']
    zc = ds['zc']
    tke = ds['tke'][t0:t1+1, :, xIdx]
    enstrophy = ds['enstrophy'][t0:t1+1, :, xIdx]
    th = ds['th'][t0:t1+1, :, xIdx]
    ds2 = xr.open_dataset(f"{DatasetDir}/{filename2}")
    wanothano = ds2['wanothano'][t0:t1+1, :, xIdx]
     
    return zc, enstrophy, tke, th

zc, enstrophy, tke, th = loadData(t0, t1, xIdx)

#%%

import vvmtools
nz, nt = 50, 721
dim_data_dict = {
    "time": (np.arange(nt)*np.timedelta64(2,'m')+np.datetime64('2024-01-01 05:00:00')).astype('datetime64[s]'),
    "height": np.arange(nz)*0.04
}
data_dict = {
    "th": th.to_numpy(),
    "enstrophy": enstrophy.to_numpy(),
    "tke": tke.to_numpy()
}
var_dims_dict = {
    "th": ("time", "height"),
    "enstrophy": ("time", "height"),
    "tke": ("time", "height")
}
attributes = {
    "th": {"units": "K", "description": "x-y mean potential temperature (t,z)"},
    "enstrophy": {"units": "1/(s^2)", "description": "x-y mean enstrophy (t,z)"},
    "tke": {"units": "(m^2)/(s^2)", "description": "x-y mean turbulent kinetic energy (t,z)"},
    "time": {"description": "Local Time"},  # Removed 'units' for time
    "height": {"units": "m", "description": "Height in grid center"},
}

# Example of creating a NetCDF file with flexible dimensions
vvmtools.analyze.create_nc_output(f"DomainProfiles-{expname}.nc", dim_data_dict, data_dict, var_dims_dict, attributes)
