import vvmtools
import numpy as np
import matplotlib.pyplot as plt
VVMTools = vvmtools.VVMTools
VVMObjects = VVMTools("/data/yhc2080/VVM/DATA/pbl_half_PU_1")

domain = (0, 1, None, None, None, None)
u = VVMObjects.get_var("u",100, domain_range=domain)
print(u)
time_list = np.arange(0,10)
u_full = VVMObjects.get_var_parallel("u",time_steps=time_list, domain_range=domain)
print(u_full.shape)

u_mean = VVMObjects.get_var_parallel("u",time_steps=time_list, domain_range=domain, cores=5, compute_mean=True, axis=(1,2))
print(u_mean)

"""
class newVVMTools(VVMTools):
    def __init__(self, case_path):
        super().__init__(case_path)

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
