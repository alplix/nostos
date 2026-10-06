import numpy as np
from galpy.potential import HernquistPotential, MiyamotoNagaiPotential, NFWPotential
from astropy import units as u

NOSTOS_G = 4.30091e-6 / (977.7922216807892**2)  # kpc^3/Msun/Myr^2
KMS_PER_KPCMYR = 977.7922216807892

bulge_GM = NOSTOS_G * 1.0e10; bulge_c = 0.7
disk_GM  = NOSTOS_G * 6.8e10; disk_a = 3.0; disk_b = 0.28
halo_GM  = NOSTOS_G * 5.4e11; halo_rs = 16.0

def nostos_bulge(R, z):
    r = np.sqrt(R*R+z*z)
    denom = r + bulge_c
    f = -bulge_GM/(denom*denom*(r if r>0 else 1e-12))
    return f*R, f*z

def nostos_disk(R, z):
    zeta = np.sqrt(z*z+disk_b*disk_b)
    azeta = disk_a+zeta
    D2 = R*R+azeta*azeta
    D3 = D2*np.sqrt(D2)
    aR = -disk_GM*R/D3
    az = -disk_GM*z*azeta/(zeta if zeta>0 else 1e-12)/D3
    return aR, az

def nostos_halo(R, z):
    r = np.sqrt(R*R+z*z)
    x = (r if r>0 else 1e-12)/halo_rs
    menc = np.log(1+x) - x/(1+x)
    r_safe = r if r>0 else 1e-12
    f = -halo_GM*menc/(r_safe**3)
    return f*R, f*z

bulge_g = HernquistPotential(amp=2*1.0e10*u.Msun, a=0.7*u.kpc)
disk_g = MiyamotoNagaiPotential(amp=6.8e10*u.Msun, a=3.0*u.kpc, b=0.28*u.kpc)
halo_g = NFWPotential(amp=5.4e11*u.Msun, a=16.0*u.kpc)

points = [(1.0,0.0),(4.0,0.0),(8.122,0.0),(8.122,0.5),(20.0,0.0)]
for R,z in points:
    print("=== R=%s, z=%s ===" % (R,z))
    for name, nf, gp in [('bulge',nostos_bulge,bulge_g), ('disk',nostos_disk,disk_g), ('halo',nostos_halo,halo_g)]:
        aR_n, az_n = nf(R,z)
        aR_g = gp.Rforce(R*u.kpc, z*u.kpc, use_physical=True) / KMS_PER_KPCMYR
        az_g = gp.zforce(R*u.kpc, z*u.kpc, use_physical=True) / KMS_PER_KPCMYR
        rel_err_R = abs(aR_n-aR_g)/abs(aR_g) if aR_g!=0 else abs(aR_n)
        rel_err_z = abs(az_n-az_g)/abs(az_g) if az_g!=0 else abs(az_n)
        print("  %-6s aR: nostos=%.6e galpy=%.6e relerr=%.2e | az: nostos=%.6e galpy=%.6e relerr=%.2e" % (name,aR_n,aR_g,rel_err_R,az_n,az_g,rel_err_z))
