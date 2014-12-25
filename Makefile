CC=g++
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
LFLAGS=-std=c++11
else
LFLAGS=-lstdc++
endif
CFLAGS=-c -Wall -g
OBJ=obj
DFLAGS=
ifdef DEBUG
	DFLAGS+=-DDEBUG=1 
endif
ifdef PARALLEL
	DFLAGS+= -fopenmp -DPARALLEL=1
endif

all: raytrace

raytrace: $(OBJ)/main.o $(OBJ)/parse.o $(OBJ)/raytrace.o $(OBJ)/image.o $(OBJ)/pixel.o $(OBJ)/EasyBMP.o
	$(CC) -g $(DFLAGS) $(LFLAGS) -o $@ $^

$(OBJ)/main.o: src/main.cpp
	$(CC) $(CFLAGS) $(DFLAGS) src/main.cpp -o $@

$(OBJ)/parse.o: src/parse.cpp
	$(CC) $(CFLAGS) $(DFLAGS) src/parse.cpp -o $@

$(OBJ)/raytrace.o: src/raytrace.cpp
	$(CC) $(CFLAGS) $(DFLAGS) src/raytrace.cpp -o $@

$(OBJ)/image.o: src/image.cpp
	$(CC) $(CFLAGS) src/image.cpp -o $@

$(OBJ)/pixel.o: src/pixel.cpp
	$(CC) $(CFLAGS) src/pixel.cpp -o $@

$(OBJ)/EasyBMP.o: src/EasyBMP.cpp
	$(CC) $(CFLAGS) src/EasyBMP.cpp -o $@

cleanObj:
	rm -rf $(OBJ)/*.o

clean:
	rm -rf $(OBJ)/*.o raytrace