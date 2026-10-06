import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic

# Real Gaia DR2 astrometry for HIP 3757 (source_id 2525688198820543360), Bailer-Jones et al.
# 2018's closest 'Oumuamua-encounter candidate -- same star tests/test_pipeline_oumuamua_hip3757.c
# and tests/test_coords_astropy_crosscheck.c use. Produces the reference numbers hardcoded in the
# latter test.
c = SkyCoord(
    ra=12.05634594268451 * u.deg, dec=-5.1359892365165924 * u.deg,
    distance=(1000.0 / 40.20064272792413) * u.pc,
    pm_ra_cosdec=182.16070593050424 * u.mas / u.yr,
    pm_dec=-120.15964959812227 * u.mas / u.yr,
    radial_velocity=26.880687935465836 * u.km / u.s,
    frame="icrs",
)
g = c.transform_to(Galactic())
cart = g.cartesian
vel = g.velocity
print("Galactic Cartesian position (pc):", cart.x.to(u.pc), cart.y.to(u.pc), cart.z.to(u.pc))
print("Galactic Cartesian velocity (km/s):", vel.d_x.to(u.km / u.s), vel.d_y.to(u.km / u.s), vel.d_z.to(u.km / u.s))
