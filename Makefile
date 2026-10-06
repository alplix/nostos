CC = gcc
CFLAGS = -O2 -std=c99 -Wall -Wextra -lm

.PHONY: test clean

test: build/test_kepler build/test_milkyway_galpy build/test_milkyway_allen_santillan build/test_coords build/test_encounter build/test_orbital_elements build/test_horizons_crosscheck build/test_montecarlo build/test_mc_search
	./build/test_kepler
	./build/test_milkyway_galpy
	./build/test_milkyway_allen_santillan
	./build/test_coords
	./build/test_encounter
	./build/test_orbital_elements
	./build/test_horizons_crosscheck
	./build/test_montecarlo
	./build/test_mc_search

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

build/test_encounter: tests/test_encounter.c src/encounter.c src/encounter.h src/potential.c src/potential.h src/integrator.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_encounter.c src/encounter.c src/potential.c -o build/test_encounter -lm

build/test_orbital_elements: tests/test_orbital_elements.c src/orbital_elements.c src/orbital_elements.h src/coords.c src/coords.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_orbital_elements.c src/orbital_elements.c src/coords.c -o build/test_orbital_elements -lm

build/test_horizons_crosscheck: tests/test_horizons_crosscheck.c src/orbital_elements.c src/orbital_elements.h src/coords.c src/coords.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_horizons_crosscheck.c src/orbital_elements.c src/coords.c -o build/test_horizons_crosscheck -lm

build/test_montecarlo: tests/test_montecarlo.c src/montecarlo.c src/montecarlo.h | build
	$(CC) $(CFLAGS) tests/test_montecarlo.c src/montecarlo.c -o build/test_montecarlo -lm

build/test_mc_search: tests/test_mc_search.c src/mc_search.c src/mc_search.h src/montecarlo.c src/encounter.c src/potential.c src/coords.c | build
	$(CC) $(CFLAGS) tests/test_mc_search.c src/mc_search.c src/montecarlo.c src/encounter.c src/potential.c src/coords.c -o build/test_mc_search -lm

clean:
	rm -rf build
