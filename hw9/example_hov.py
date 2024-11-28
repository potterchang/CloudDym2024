import numpy as np
from vvmtools.plot import DataPlotter
import matplotlib.pyplot as plt
import xarray as xr

# prepare expname and data coordinate
expname  = 'pbl_PU_s1'
nx = 128; x = np.arange(nx)*0.2
ny = 128; y = np.arange(ny)*0.2
nz = 50;  z = np.arange(nz)*0.04
nt = 721; t = np.arange(nt)*np.timedelta64(2,'m')+np.datetime64('2024-01-01 05:00:00')

#Load Data
ds = xr.open_dataset(f"/data/yhc2080/CloudDym2024/DATA/Pollutants_{expname}-000000-000720.nc")
data_xt2d = ds['NO2'].isel(z=1) +  ds['NO'].isel(z=1)


#%%
# create dataPlotter class
figpath           = './fig/'
data_domain       = {'x':x, 'y':y, 'z':z, 't':t}
data_domain_units = {'x':'km', 'y':'km', 'z':'km', 't':'LocalTime'}
dplot = DataPlotter(expname, figpath, data_domain, data_domain_units)

# TODO: Change to your own data
np.random.seed(0)
#data_xt2d  = np.random.normal(0, 0.1, size=(nt,nx))

fig, ax, cax = dplot.draw_xt(data = data_xt2d,
                                levels = np.arange(0,30.001,2),
                                extend = 'max',
                                cmap_name = 'OrRd',
                                title_left  = 'NOx (ppb)',
                                title_right = f'Domain',
                                figname     = f'NOx_{expname}.png',
                               )
plt.show()
