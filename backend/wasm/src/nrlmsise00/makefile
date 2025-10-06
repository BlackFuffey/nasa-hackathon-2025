CFLAGS = -Wall -g -DINLINE
MFLAGS = -lm
WASM = emcc $(CFLAGS)

nrlmsise00.wasm :	nrlmsise-00.c nrlmsise-00.h
	$(WASM) nrlmsise-00.h -o nrlmsise00.wasm

clean   :
	rm -rf  nrlmsise-test nrlmsise-00.o nrlmsise-00_test.o \
                   nrlmsise-00_data.o
