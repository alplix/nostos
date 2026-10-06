CC = gcc
CFLAGS = -O2 -std=c99 -Wall -Wextra -lm

.PHONY: test clean

test: build/test_kepler build/test_milkyway build/test_coords
	./build/test_kepler
	./build/test_milkyway
	./build/test_coords

build:
	mkdir -p build

build/test_kepler: tests/test_kepler.c src/potential.c src/potential.h src/integrator.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_kepler.c src/potential.c -o build/test_kepler -lm

build/test_milkyway: tests/test_milkyway.c src/potential.c src/potential.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_milkyway.c src/potential.c -o build/test_milkyway -lm

build/test_coords: tests/test_coords.c src/coords.c src/coords.h src/vec3.h src/units.h | build
	$(CC) $(CFLAGS) tests/test_coords.c src/coords.c -o build/test_coords -lm

clean:
	rm -rf build
