import vvmtools
import os
import re
import traceback
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
VVMTools = vvmtools.VVMTools
VVMObjects = VVMTools("/data/yhc2080/VVM/DATA/pbl_half_PU_1")

"""
domain = (0, 1, None, None, None, None)
u = VVMObjects.get_var("u",100, domain_range=domain)
print(u)
time_list = np.arange(0,10)
u_full = VVMObjects.get_var_parallel("u",time_steps=time_list, domain_range=domain)
print(u_full.shape)

u_mean = VVMObjects.get_var_parallel("u",time_steps=time_list, domain_range=domain, cores=5, compute_mean=True, axis=(1,2))
print(u_mean)
"""

#%%
class myVVMTools(VVMTools):
    def __init__(self, case_path):
        super().__init__(case_path)
    
    def get_var_file_path(self,var,time=0):
        
        variable_type = self.get_variable_file_type(var)
        if self.DEBUGMODE:
            print(f"Variable type: {variable_type}")
        if variable_type == "Variable not found":
            print(f"Variable {var} not found in the dataset.")
            return None       

        if variable_type == "TOPO":
            # Special case for TOPO variables, always in TOPO.nc
            topo_file = os.path.join(self.CASEPATH, 'TOPO.nc')
            file_path = topo_file
        else:            
            # Construct the expected filename pattern for the given variable and time
            file_pattern = f"{variable_type}-{'{:06d}'.format(time)}.nc"
            regex_pattern = f".*{file_pattern}$"  # Convert glob pattern to regex
            if self.DEBUGMODE:
                print(f"Regex Pattern: {regex_pattern}")
            # Search for the file in the case path
            for root, dirs, files in os.walk(self.CASEPATH):
                for filename in files:
                    if re.match(regex_pattern, filename):
                        if self.DEBUGMODE:
                            print(f"File found: {filename}")
                        # Uncomment the following block to open and read the file
                        file_path = os.path.join(root, filename)
        return file_path
    
       
    def get_var_regrid(self, 
                var, 
                time, 
                domain_range=(None, None, None, None, None, None), # (k1, k2, j1, j2, i1, i2)
                TOPOmask=False, fill_value=None, TOPOmethod="nearest",
                numpy=False, 
                compute_mean=False, 
                axis=None):   
        nz, ny, nx = self.DIM['zc'].size, self.DIM['yc'].size, self.DIM['xc'].size
        domain_full_range = [0, nz, 0, ny, 0, nx]
        if type(domain_range)==tuple:
            domain_range = list(domain_range)
        for h in range(6):
            if domain_range[h] is None:
                domain_range[h] = domain_full_range[h]
        k1, k2, j1, j2, i1, i2 = domain_range
        ResultSlicerList = [slice(0,1),slice(k1, k2),slice(j1,j2),slice(i1,i2)]
        try:
            filepath = mytools.get_var_file_path("u",time)
            ds = xr.open_dataset(filepath)
            if TOPOmask == True:
                topopath = mytools.get_var_file_path("topo",0)
                dstopo = xr.open_dataset(topopath)
                try:
                    mask_full = dstopo['mask'].astype(bool)
                except: # create a new mask
                    TOPOindex = dstopo['topo'] -1 
                    TOPOindex[TOPOindex==-1] = 0
                    levelindex = np.arange(0,nz)
                    TOPO3D = np.expand_dims(TOPOindex, axis=0)
                    level3D = np.expand_dims(levelindex, axis=(1,2))
                    TMask = level3D <= TOPO3D # True for terrain; False for air
                    mask_full = ~TMask
                # TOPO mask value setting
                mask_values = {'nearest': np.nan, 'linear0': 0}
                if fill_value is None:
                    fill_value = mask_values[TOPOmethod]
            
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

            
            # Regridding
            if var == "eta_2":
                var = "eta"
            varName = var
            Var_ori_sliced = ds[varName][tuple(ResultSlicerList)]
            varshape = ds[varName].shape
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
                Var_rolled_sliced = ds[varName][tuple(RolledSlicerList)]
                # TOPO mask
                if TOPOmask == True:
                    TOPOmask_ori = ~mask_full
                    TOPOmask_rolled = np.roll(~mask_full,mask_roll[varName],axis=axis_roll[varName])
                    TOPOmask_ori_sliced = np.expand_dims(TOPOmask_ori[tuple(ResultSlicerList[-3:])],axis=0)
                    TOPOmask_rolled_sliced = np.expand_dims(TOPOmask_rolled[tuple(ResultSlicerList[-3:])],axis=0)
                    if np.isscalar(Var_ori_sliced):
                        if TOPOmask_ori_sliced:
                            Var_ori_sliced = fill_value
                    else:
                        Var_ori_sliced= Var_ori_sliced.where(~TOPOmask_ori_sliced, fill_value)
                    if np.isscalar(Var_rolled_sliced):
                        if TOPOmask_rolled_sliced:
                            Var_rolled_sliced = fill_value
                    else:
                        Var_rolled_sliced[TOPOmask_rolled_sliced] = fill_value
    
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
                Var_rolled_sliced1 = ds[varName][tuple(RolledSlicerList1)]
                Var_rolled_sliced2 = ds[varName][tuple(RolledSlicerList2)]
                Var_rolled_sliced12 = ds[varName][tuple(RolledSlicerList12)]
                if TOPOmask == True:
                    TOPOmask_ori = ~mask_full
                    TOPOmask_rolled1 = np.roll(~mask_full,mask_roll1[varName],axis=axis_roll1[varName])
                    TOPOmask_rolled2 = np.roll(~mask_full,mask_roll2[varName],axis=axis_roll2[varName])
                    TOPOmask_rolled12 = np.roll(np.roll(~mask_full,mask_roll1[varName],axis=axis_roll1[varName]),mask_roll2[varName],axis=axis_roll2[varName])   
                    Var_original[TOPOmask_ori] = fill_value
                    Var_rolled1[TOPOmask_rolled1] = fill_value
                    Var_rolled2[TOPOmask_rolled2] = fill_value
                    Var_rolled12[TOPOmask_rolled12] = fill_value
    
                # Calculate regridded value
                OUT = np.nanmean([Var_ori_sliced,Var_rolled_sliced1,Var_rolled_sliced2,Var_rolled_sliced12],axis=0)    
                
            else:
                if TOPOmask == True:
                    Var_ori_sliced[~mask] = fill_value
                OUT = Var_ori_sliced
            ds.close()
            if numpy:
                data = np.squeeze(OUT.to_numpy())

                if compute_mean and axis is not None:
                    return np.mean(data, axis=axis)
                elif compute_mean:
                    return np.mean(data)
                else:
                    return data
    
            return OUT
    
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            traceback.print_exc()
            return None  

    

    def cal_TKE(self, t):
        u = np.squeeze(self.get_var("u",t, numpy=True))
        v = np.squeeze(self.get_var("v",t, numpy=True))
        w = np.squeeze(self.get_var("w",t, numpy=True))
        return np.nanmean(1/2*(u**2+v**2+w**2), axis=(1,2))

    def cal_Enstrophy(self, t):
        if self.VARTYPE["eta"]=="Dynamic":
            eta = np.squeeze(self.get_var("eta"),t,numpy=True)
        else:
            eta = np.squeeze(self.get_var("eta_2",t, numpy=True))
        zeta = np.squeeze(self.get_var("zeta",t, numpy=True))
        xi = np.squeeze(self.get_var("xi",t, numpy=True))
        return np.nanmean((eta**2+zeta**2+xi**2), axis=(1,2))
        
mytools = myVVMTools("/data/yhc2080/VVM/DATA/pbl_half_PU_1")

filepath = mytools.get_var_file_path("u",100)
#%%
ds = xr.open_dataset(filepath)

#%%
newtools = newVVMTools("/data/yhc2080/VVM/DATA/pbl_half_PU_1")
TKE_all = newtools.func_time_parallel(newtools.cal_TKE, np.arange(0,10),cores=5)

Ens_all = newtools.func_time_parallel(newtools.cal_Enstrophy, np.arange(0,10),cores=5)

print(Ens_all.shape)

t = np.arange(10)
z = newtools.DIM["zc"]
tt, zz = np.meshgrid(t,z)
plt.contourf(t,z,Ens_all.T)
plt.show()
"""
