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

#%%

def find_BL_boundary(var, howToSearch, threshold=0.01, zc=z):
    nt, nz = var.shape[0], var.shape[1]

    if howToSearch == "th_plus05K":
        # search the index where the var is closest to (first reach) var_sfc+0.5K
        thbot = var[:,1] 
        P05K = thbot + 0.5
        ind_posi = np.array((var - P05K) >= 0)
        adj = np.sum(ind_posi, axis=1) == 0
        ind_BLH1 = np.argmax(ind_posi, axis=1) - adj
        ind_BLH1 = np.ma.array(ind_BLH1, mask=ind_BLH1<0)
        BLH1 = np.take(zc, ind_BLH1)
        return BLH1
    elif howToSearch == "dthdz":
        # search the index where the dvar/dz is the largest
        dthdz = np.gradient(var, zc, axis=1)
        ind_BLH2 = np.argmax(dthdz, axis=1)
        BLH2 = np.take(zc, ind_BLH2)
        return BLH2
    elif howToSearch == "threshold":
        # search the index where var is first change from exceed over to less than the threshold value
        ind_posi = np.array((var - threshold) >= 0)
        BLH3 = np.zeros(nt)
        for t in range(nt):
            arr = ind_posi[t]
            ind_BLH3 = 0
            for i in range(1, len(arr)):
                if not arr[i] and arr[i - 1]:
                    ind_BLH3 = i - 1
                    continue
            BLH3[t] = np.take(zc, ind_BLH3)
        
        return BLH3
    elif howToSearch == "wth":
        BLHs = np.zeros((3,nt))
        # 3 wth definition: 
        # first boundary from + to -
        ind_posi = np.array((var) < 0)
        adj = np.sum(ind_posi, axis=1) == 0
        ind_BLHw1 = np.argmax(ind_posi, axis=1) - adj
        ind_BLHw1 = np.ma.array(ind_BLHw1, mask=ind_BLHw1<0)
        BLHs[0, :] = np.take(zc, ind_BLHw1)
        # negative most
        BLHs[1, :] =  np.take(zc, np.argmin(var, axis=1))
        # first boundary from - to +
        for t in range(nt):
            ind_posi = np.array((var[t]) > 0)
            ind_posi[:np.argmin(var[t])] = False
            adj = np.sum(ind_posi) == 0
            ind_BLHw3 = np.argmax(ind_posi) - adj
            ind_BLHw3 = np.ma.array(ind_BLHw3, mask=ind_BLHw3<0)
            BLHs[2, t] = np.take(zc, ind_BLHw3)          
            
        # Exception
        for t in range(nt):                
            if np.nanmax(var[t])<1e-3: # very small wth
                BLHs[:,t] = 0
            elif np.nanmax(var[t, np.argmin(var[t]):]) < 1e-3: # no apparent boundary from - to +
                BLHs[2,t] = 0
        return BLHs[0], BLHs[1], BLHs[2]
            
    else:
        print("Without this searching approach")
        return np.nan

def find_BL_boundary2(var, howToSearch, threshold=0.01, zc=z):
    """
    Calculate the boundary layer height based on specified search criteria.

    :param var: Array of the variable used to define the boundary layer (e.g., potential temperature θ).
    :param howToSearch: String defining the boundary search method. Options are:
                        - "th_plus05K": Height where θ exceeds surface value by 0.5K.
                        - "dthdz": Height where the vertical gradient of θ is maximum.
                        - "threshold": Height where `var` first exceeds the threshold value.
                        - "wth": Boundary levels based on transitions in vertical heat flux.
    :param threshold: Threshold value used in the "threshold" and "wth" methods (default is 0.01).
    :return: Array of boundary layer heights for each time step or, for "wth", an array with lower, mid, and upper boundaries.
    """
    
    var = np.array(var)
    
    if howToSearch == "th_plus05K":
        # Find the height where theta exceeds surface value by 0.5K
        th_find = var - (var[:, 0].reshape(var.shape[0], 1) + 0.5)

        # Identify boundary layer height for each time step
        h_BL_th_plus05 = []
        for t in range(len(th_find)):
            h_BL_th_plus05.append(zc[np.where(th_find[t] < threshold)[0][-1]])
        return np.array(h_BL_th_plus05)

    elif howToSearch == "dthdz":
        # Compute vertical gradient of theta (dθ/dz)
        dth_dz = (var[:, 1:] - var[:, :-1]) / (zc[1:] - zc[:-1])

        # Determine height of maximum gradient for each time step
        #zc = zc[1:]
        h_BL_dthdz = []
        for t in range(len(dth_dz)):
            h_BL_dthdz.append(zc[np.argmax(dth_dz[t])])
        return np.array(h_BL_dthdz)
    elif howToSearch == "threshold":
        # Identify height where 'var' crosses the threshold
        positive_mask = np.swapaxes(var, 0, 1) > threshold
        index = np.argwhere(positive_mask)

        # Map each time step to the highest index where the threshold is crossed
        k, i = index[:, 0], index[:, 1]
        i_k_map = {}

        for ii in range(len(i)):
            if i[ii] not in i_k_map:
                i_k_map[i[ii]] = k[ii]
            else:
                if k[ii] > i_k_map[i[ii]]:
                    i_k_map[i[ii]] = k[ii]

        # Compute boundary layer height for each time step based on threshold crossings
        #zc = zc[1:]
        h = np.zeros(len(var))
        for key, value in i_k_map.items():
            h[key] = zc[value]
        return np.array(h)

    elif howToSearch == "wth":
        # Identify lower, middle, and upper boundaries based on vertical heat flux transitions
        mask = np.where(var>=0,1,-1)
        zc_wth = np.zeros((3,var.shape[0]))
        #zc = zc[1:]
        for t in range(var.shape[0]):
            temp = mask[t,1:] - mask[t,:-1]

            # Default zero boundary if w'θ' values are low
            if np.max(var[t])<threshold:
                k_lower = k_mid = k_upper = 0

            # Find lower, mid, and upper boundaries based on transitions
            else:
                # Find lower boundary
                try:
                    k_lower = np.argwhere(temp==-2)[0][0]
                except:
                    k_lower = 0

                # Find mid boundary
                k_mid = np.argmin(var[t])

                # Find upper boundary
                try:
                    k_upper = np.argwhere(temp==2)[0][0]
                except:
                    k_upper = 0

            # Zero out upper boundary if max value beyond mid-boundary is low
            if np.max(var[t, np.argmin(var[t]):]) < threshold:
                k_upper = 0

            # Store the boundaries for this time step
            zc_wth[0,t], zc_wth[1,t], zc_wth[2,t] = zc[k_lower], zc[k_mid], zc[k_upper]
        return zc_wth

    else:
        print("Without this searching approach")


# Read regional averaged data
DatasetDir = "../DATA"
filename = f"Turbulence_{expname}-000000-000720.nc"
filename2 = f"EddyFlux_{expname}-000000-000720.nc"

t0 = 0
t1 = 721
xIdx= 2
trIdx = 2
loc_list = ["pasture","urban", "domain"]
slice_list = [slice(0,64), slice(64,128), slice(0,128)]
islice = slice_list[xIdx]
tr_list = ["tracer_sfc", "tracer_750", "tracer_1500"]

# Load Tracer Data
ds = xr.open_dataset(f"/data/yhc2080/CloudDym2024/DATA/Pollutants_{expname}-000000-000720.nc")
data_zt2d = ds[f'tr{str(trIdx+1).zfill(2)}'].isel(x=islice).mean(dim="x").T
data_zt2d = data_zt2d / np.max(data_zt2d)

def loadData(t0, t1, xIdx):
    ds = xr.open_dataset(f"{DatasetDir}/{filename}")
    time = ds['time']
    zc = ds['zc']
    tke = ds['tke'][t0:t1+1, :, xIdx]
    enstrophy = ds['enstrophy'][t0:t1+1, :, xIdx]
    th = ds['th'][t0:t1+1, :, xIdx]
    ds2 = xr.open_dataset(f"{DatasetDir}/{filename2}")
    wanothano = ds2['wanothano'][t0:t1+1, :, xIdx]
    
    wanothano[:,0] = wanothano[:,1]
    
    BLH1 = find_BL_boundary(th, "th_plus05K")
    BLH2 = find_BL_boundary(th, "dthdz")
    BLH3 = find_BL_boundary(tke, "threshold", threshold=0.08)
    BLH4 = find_BL_boundary(enstrophy, "threshold", threshold=1e-5)
    BLH5, BLH6, BLH7 = find_BL_boundary(np.array(wanothano), "wth")

    BLHs = [BLH1, BLH2, BLH3, BLH4, BLH5, BLH6, BLH7]

    return wanothano, BLHs, zc, enstrophy, tke

wanothano, BLHs, zc, enstrophy, tke = loadData(t0, t1, xIdx)
BLHtext = [r"$\theta_{sfc}+0.5K$",r"max $d\theta/dz$", r"TKE=0.08","Enstrophy=1e-5",
           r"top($\overline{w'\theta'}+$)",r"min($\overline{w'\theta'}$)",r"top($\overline{w'\theta'}-$)"]
pblh_dicts = {f"{BLHtext[i]}": BLHs[i] for i in range(len(BLHs))}


# create dataPlotter class
figpath           = './fig/'
data_domain       = {'x':x, 'y':y, 'z':z, 't':t}
data_domain_units = {'x':'km', 'y':'km', 'z':'km', 't':'LocalTime'}
dplot = DataPlotter(expname, figpath, data_domain, data_domain_units)

# TODO: Use your tracer data, noted that you need to normalize the data by its maximum
print(trIdx, xIdx)
fig, ax, cax = dplot.draw_zt(data = data_zt2d,
                            levels = np.arange(0,1.1,0.1),
                            extend = 'neither',
                            pblh_dicts=pblh_dicts,
                            title_left  = f'{tr_list[trIdx]} (normalized)',
                            title_right = f'{loc_list[xIdx]}',
                            cmap_name   = 'Greys',
                            figname     = 'test_pbl.png',
                      )


### If you want to delete the legend, turn this block on
ax.get_legend().remove()
plt.savefig(f'{figpath}/{tr_list[trIdx]}_{expname}_{loc_list[xIdx]}.png', dpi=200)
plt.close('all')
###

plt.show()


#%%

fig, ax, cax = dplot.draw_zt(data = wanothano.T,
                            levels = np.arange(-0.04,0.04,0.01),
                            extend = 'both',
                            pblh_dicts=pblh_dicts,
                            title_left  = "w'th'", # f'{tr_list[trIdx]} (normalized)',
                            title_right = f'{loc_list[xIdx]}',
                            cmap_name   = 'coolwarm',
                            figname     = "w'th'2.png",
                      )
