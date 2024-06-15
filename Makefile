# Compiler flags
CFLAGS = -O3 -fPIC -shared
LDFLAGS = -lm -lfftw3f -fopenmp -lfftw3f_omp

# Rule to build the shared library
shared_lib: Reion.so

Reion.so: Reionfuncs.c
	gcc $(CFLAGS) Reionfuncs.c -o Reion.so $(LDFLAGS)

# Rule to clean the build files
clean:
	rm -f Reion.so

.PHONY: shared_lib clean

