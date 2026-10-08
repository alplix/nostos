CC = gcc
PYTHON ?= python
CFLAGS = -O2 -std=c99 -Wall -Wextra -lm

.PHONY: test clean tools

test: build/test_kepler build/test_milkyway_galpy build/test_milkyway_allen_santillan build/test_coords build/test_coords_astropy_crosscheck build/test_encounter build/test_orbital_elements build/test_horizons_crosscheck build/test_montecarlo build/test_mc_search build/test_pipeline_oumuamua_hip3757 build/test_pipeline_borisov_ross573 build/test_scan_logic build/test_nostos_io build/test_star_mc build/test_parallax_zeropoint build/test_csv_range
	./build/test_kepler
	./build/test_milkyway_galpy
	./build/test_milkyway_allen_santillan
	./build/test_coords
	./build/test_coords_astropy_crosscheck
	./build/test_encounter
	./build/test_orbital_elements
	./build/test_horizons_crosscheck
	./build/test_montecarlo
	./build/test_mc_search
	./build/test_pipeline_oumuamua_hip3757
	./build/test_pipeline_borisov_ross573
	./build/test_scan_logic
	./build/test_nostos_io
	./build/test_star_mc
	./build/test_parallax_zeropoint
	./build/test_csv_range
	$(PYTHON) tests/test_make_work_units.py

build:
	mkdir -p build

build/test_kepler: tests/test_kepler.c src/potential.c src/potential.h src/integrator.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_kepler.c src/potential.c -o build/test_kepler -lm

build/test_milkyway_galpy: tests/test_milkyway_galpy.c src/potential.c src/potential.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_milkyway_galpy.c src/potential.c -o build/test_milkyway_galpy -lm

build/test_milkyway_allen_santillan: tests/test_milkyway_allen_santillan.c src/potential.c src/potential.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_milkyway_allen_santillan.c src/potential.c -o build/test_milkyway_allen_santillan -lm

build/test_coords: tests/test_coords.c src/coords.c src/coords.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_coords.c src/coords.c -o build/test_coords -lm

build/test_coords_astropy_crosscheck: tests/test_coords_astropy_crosscheck.c src/coords.c src/coords.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_coords_astropy_crosscheck.c src/coords.c -o build/test_coords_astropy_crosscheck -lm

build/test_encounter: tests/test_encounter.c src/encounter.c src/encounter.h src/potential.c src/potential.h src/integrator.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_encounter.c src/encounter.c src/potential.c -o build/test_encounter -lm

build/test_orbital_elements: tests/test_orbital_elements.c src/iso_solutions.h src/orbital_elements.c src/orbital_elements.h src/coords.c src/coords.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_orbital_elements.c src/orbital_elements.c src/coords.c -o build/test_orbital_elements -lm

build/test_horizons_crosscheck: tests/test_horizons_crosscheck.c src/orbital_elements.c src/orbital_elements.h src/coords.c src/coords.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_horizons_crosscheck.c src/orbital_elements.c src/coords.c -o build/test_horizons_crosscheck -lm

build/test_montecarlo: tests/test_montecarlo.c src/montecarlo.c src/montecarlo.h | build
	$(CC) $(CFLAGS) tests/test_montecarlo.c src/montecarlo.c -o build/test_montecarlo -lm

build/test_mc_search: tests/test_mc_search.c src/mc_search.c src/mc_search.h src/montecarlo.c src/encounter.c src/potential.c src/coords.c | build
	$(CC) $(CFLAGS) tests/test_mc_search.c src/mc_search.c src/montecarlo.c src/encounter.c src/potential.c src/coords.c -o build/test_mc_search -lm

build/test_pipeline_oumuamua_hip3757: tests/test_pipeline_oumuamua_hip3757.c src/iso_solutions.h src/orbital_elements.h src/orbital_elements.c src/coords.c src/encounter.c src/potential.c | build
	$(CC) $(CFLAGS) tests/test_pipeline_oumuamua_hip3757.c src/orbital_elements.c src/coords.c src/encounter.c src/potential.c -o build/test_pipeline_oumuamua_hip3757 -lm

build/test_pipeline_borisov_ross573: tests/test_pipeline_borisov_ross573.c src/iso_solutions.h src/orbital_elements.h src/orbital_elements.c src/coords.c src/encounter.c src/potential.c | build
	$(CC) $(CFLAGS) tests/test_pipeline_borisov_ross573.c src/orbital_elements.c src/coords.c src/encounter.c src/potential.c -o build/test_pipeline_borisov_ross573 -lm

build/test_scan_logic: tests/test_scan_logic.c src/scan_logic.c src/scan_logic.h src/encounter.c src/potential.c | build
	$(CC) $(CFLAGS) tests/test_scan_logic.c src/scan_logic.c src/encounter.c src/potential.c -o build/test_scan_logic -lm

build/test_nostos_io: tests/test_nostos_io.c src/nostos_io.c src/nostos_io.h src/coords.c | build
	$(CC) $(CFLAGS) tests/test_nostos_io.c src/nostos_io.c src/coords.c -o build/test_nostos_io -lm

build/test_star_mc: tests/test_star_mc.c src/star_mc.c src/star_mc.h src/montecarlo.c src/montecarlo.h | build
	$(CC) $(CFLAGS) tests/test_star_mc.c src/star_mc.c src/montecarlo.c -o build/test_star_mc -lm

tools: build/gaia_scan build/oumuamua_velocity_covariance build/mc_vs_bestfit_ross573 build/mc_shortlist

build/gaia_scan: src/gaia_scan.c src/iso_solutions.h src/orbital_elements.h src/coords.c src/orbital_elements.c src/encounter.c src/potential.c src/scan_logic.c src/parallax_zeropoint.c src/parallax_zeropoint.h src/csv_range.c src/csv_range.h | build
	$(CC) $(CFLAGS) -Isrc src/gaia_scan.c src/coords.c src/orbital_elements.c src/encounter.c src/potential.c src/scan_logic.c src/parallax_zeropoint.c src/csv_range.c -o build/gaia_scan -lm

build/oumuamua_velocity_covariance: tools/oumuamua_velocity_covariance.c src/iso_solutions.h src/orbital_elements.c src/coords.c | build
	$(CC) $(CFLAGS) tools/oumuamua_velocity_covariance.c src/orbital_elements.c src/coords.c -o build/oumuamua_velocity_covariance -lm

build/mc_vs_bestfit_ross573: tools/mc_vs_bestfit_ross573.c src/iso_solutions.h src/star_mc.c src/star_mc.h src/orbital_elements.c src/coords.c src/encounter.c src/potential.c src/montecarlo.c | build
	$(CC) $(CFLAGS) tools/mc_vs_bestfit_ross573.c src/star_mc.c src/orbital_elements.c src/coords.c src/encounter.c src/potential.c src/montecarlo.c -o build/mc_vs_bestfit_ross573 -lm

build/mc_shortlist: tools/mc_shortlist.c src/star_mc.c src/star_mc.h src/iso_solutions.h src/orbital_elements.c src/coords.c src/encounter.c src/potential.c src/montecarlo.c src/scan_logic.c src/parallax_zeropoint.c src/parallax_zeropoint.h src/csv_range.c src/csv_range.h | build
	$(CC) $(CFLAGS) -Isrc tools/mc_shortlist.c src/star_mc.c src/orbital_elements.c src/coords.c src/encounter.c src/potential.c src/montecarlo.c src/scan_logic.c src/parallax_zeropoint.c -o build/mc_shortlist -lm

build/test_parallax_zeropoint: tests/test_parallax_zeropoint.c src/parallax_zeropoint.c src/parallax_zeropoint.h | build
	$(CC) $(CFLAGS) tests/test_parallax_zeropoint.c src/parallax_zeropoint.c -o build/test_parallax_zeropoint -lm

build/test_csv_range: tests/test_csv_range.c src/csv_range.c src/csv_range.h | build
	$(CC) $(CFLAGS) tests/test_csv_range.c src/csv_range.c -o build/test_csv_range -lm

clean:
	rm -rf build
