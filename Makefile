# Compiler flags
CFLAGS = -O3 -fPIC 
LDFLAGS = -lm -lfftw3f

# Rule to build the shared library
shared_lib: Reion.so

# Rule to build the object file
s.o: Reionfuncs.c
	gcc $(CFLAGS) -c Reionfuncs.c -o s.o

# Rule to link the shared library
Reion.so: s.o
	gcc -shared $(CFLAGS) -o Reion.so s.o $(LDFLAGS)

# Rule to clean the build files
clean:
	rm -f s.o Reion.so

.PHONY: shared_lib clean

