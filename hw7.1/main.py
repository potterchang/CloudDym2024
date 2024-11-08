import numpy as np
import xarray as xr
from plottools import dataPlotters
import matplotlib.pyplot as plt

# TODO 1: import your vvmtools
from vvmtools import VVMTools

# TODO 2: change the expname to your experiment name
# prepare expname and data coordinate
expname  = 'pbl_half_PU_uarea_2'
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
loc_list = ["Pasture","Urban", "Full Domain"]

vtls = VVMTools(case_path=f"/data/yhc2080/VVM/DATA/{expname}")

def loadData(t0, t1, xIdx):
    ds = xr.open_dataset(f"{DatasetDir}/{filename}")
    time = ds['time']
    zc = ds['zc']
    tke = ds['tke'][t0:t1+1, :, xIdx]
    enstrophy = ds['enstrophy'][t0:t1+1, :, xIdx]
    th = ds['th'][t0:t1+1, :, xIdx]
    ds2 = xr.open_dataset(f"{DatasetDir}/{filename2}")
    wanothano = ds2['wanothano'][t0:t1+1, :, xIdx]
     
    BLH1 = vtls.find_BL_boundary(th, "th_plus05K")
    BLH2 = vtls.find_BL_boundary(th, "dthdz")
    BLH3 = vtls.find_BL_boundary(tke, "threshold", threshold=0.02)
    BLH4 = vtls.find_BL_boundary(enstrophy, "threshold", threshold=2e-5)
    BLH5, BLH6, BLH7 = vtls.find_BL_boundary(np.array(wanothano), "wth")

    BLHs = [BLH1, BLH2, BLH3, BLH4, BLH5, BLH6, BLH7]

    return wanothano, BLHs, zc

wanothano, BLHs, zc = loadData(t0, t1, xIdx)

# read or create data
data_zt2d  = wanothano.T
BLHtext = [r"$\theta_{sfc}+0.5K$",r"max $d\theta/dz$", r"TKE=0.02","Enstrophy=2e-5",
           r"top($\overline{w'\theta'}+$)",r"min($\overline{w'\theta'}$)",r"top($\overline{w'\theta'}-$)"]
pblh_dicts = {f"{BLHtext[i]}": BLHs[i]/1000 for i in range(len(BLHs))}

# TODO 4: change the figpath to your figure path
# create dataPlotter class
figpath           = './fig/'
data_domain       = {'x':x, 'y':y, 'z':z, 't':t}
data_domain_units = {'x':'km', 'y':'km', 'z':'km', 't':'LocalTime'}
dplot = dataPlotters(expname, figpath, data_domain, data_domain_units)

# draw z-t diagram
# input data dimension is (nz, nt)
# [output] figure, axis, colorbar axis

# TODO 5: change the levels to your data range, 
#         add pblh_dicts for your seven lines, 
#         change the title_left and title_right,
#         change figname for output file name.
fig, ax, cax = dplot.draw_zt(data = data_zt2d, \
                             levels = np.arange(-0.05,0.051,0.01), \
                             extend = 'both', \
                             pblh_dicts=pblh_dicts,\
                             title_left  = r"$\overline{w'\theta'}$ & PBL Height", \
                             title_right = f'{loc_list[xIdx]}', \
                             figname     = f'Wth-{expname}-{loc_list[xIdx]}.png',\
                      )

#plt.show()

