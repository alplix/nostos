"""Quick smoke test: confirm we can query the Gaia DR3 archive for full
6D phase-space sources (position + velocity, i.e. sources with a measured
radial_velocity), and see what the real column set / row count looks like.
Not the real bulk-download script -- just validating access and format."""
from astroquery.gaia import Gaia

query = """
SELECT TOP 10
    source_id, ra, dec, parallax, parallax_error,
    pmra, pmra_error, pmdec, pmdec_error,
    radial_velocity, radial_velocity_error,
    phot_g_mean_mag
FROM gaiadr3.gaia_source
WHERE radial_velocity IS NOT NULL
  AND parallax IS NOT NULL AND parallax > 0
"""

job = Gaia.launch_job_async(query)
results = job.get_results()
print(results)
print()
print("columns:", results.colnames)

# Also get a real count of the full 6D subset size.
count_query = """
SELECT COUNT(*) as n
FROM gaiadr3.gaia_source
WHERE radial_velocity IS NOT NULL
  AND parallax IS NOT NULL AND parallax > 0
"""
count_job = Gaia.launch_job_async(count_query)
count_results = count_job.get_results()
print()
print("Full 6D (RV + positive parallax) source count:", count_results['n'][0])
