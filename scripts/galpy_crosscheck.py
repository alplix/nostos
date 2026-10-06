from galpy.potential import PlummerPotential, MiyamotoNagaiPotential
from astropy import units as u

# Allen & Santillan (1991) bulge (Plummer) and disk (Miyamoto-Nagai)
# parameters -- see potential.h/potential.c's comments for sources. The halo
# has no simple galpy equivalent (a generalized truncated enclosed-mass law,
# not a stock galpy potential class) -- see tests/test_milkyway_allen_santillan.c
# for how that component is validated instead (continuity + numerical-derivative
# self-consistency, since there's no independent reference implementation to
# compare against here).
NOSTOS_G = 4.30091e-6 / (977.7922216807892**2)  # kpc^3/Msun/Myr^2
KMS_PER_KPCMYR = 977.7922216807892
MGAL_TO_MSUN = 2.325e7

bulge_GM = NOSTOS_G * (606.0 * MGAL_TO_MSUN); bulge_b = 0.3873
disk_GM  = NOSTOS_G * (3690.0 * MGAL_TO_MSUN); disk_a = 5.3178; disk_b = 0.2500

def nostos_bulge(R, z):
    import math
    r = math.sqrt(R*R+z*z)
    denom = (r*r + bulge_b*bulge_b) ** 1.5
    f = -bulge_GM/denom
    return f*R, f*z

def nostos_disk(R, z):
    import math
    zeta = math.sqrt(z*z+disk_b*disk_b)
    azeta = disk_a+zeta
    D2 = R*R+azeta*azeta
    D3 = D2*math.sqrt(D2)
    aR = -disk_GM*R/D3
    az = -disk_GM*z*azeta/(zeta if zeta>0 else 1e-12)/D3
    return aR, az

# galpy's PlummerPotential amp convention: Phi = -amp/sqrt(r^2+b^2) directly
# (unlike Hernquist's amp=2GM; confirmed by matching the plain Plummer
# formula both use).
bulge_g = PlummerPotential(amp=606.0*MGAL_TO_MSUN*u.Msun, b=0.3873*u.kpc)
disk_g = MiyamotoNagaiPotential(amp=3690.0*MGAL_TO_MSUN*u.Msun, a=5.3178*u.kpc, b=0.2500*u.kpc)

points = [(1.0,0.0),(4.0,0.0),(8.122,0.0),(8.122,0.5),(20.0,0.0)]
for R,z in points:
    print("=== R=%s, z=%s ===" % (R,z))
    for name, nf, gp in [('bulge',nostos_bulge,bulge_g), ('disk',nostos_disk,disk_g)]:
        aR_n, az_n = nf(R,z)
        aR_g = gp.Rforce(R*u.kpc, z*u.kpc, use_physical=True) / KMS_PER_KPCMYR
        az_g = gp.zforce(R*u.kpc, z*u.kpc, use_physical=True) / KMS_PER_KPCMYR
        rel_err_R = abs(aR_n-aR_g)/abs(aR_g) if aR_g!=0 else abs(aR_n)
        rel_err_z = abs(az_n-az_g)/abs(az_g) if az_g!=0 else abs(az_n)
        print("  %-6s aR: nostos=%.6e galpy=%.6e relerr=%.2e | az: nostos=%.6e galpy=%.6e relerr=%.2e" % (name,aR_n,aR_g,rel_err_R,az_n,az_g,rel_err_z))
