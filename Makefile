CC=g++
CFLAGS=-c -Wall -std=c++11
OBJ=obj
DFLAGS=
ifdef DEBUG
	DFLAGS+=-g -DDEBUG=1 
endif
ifdef PARALLEL
	DFLAGS+= -fopenmp -DPARALLEL=1
endif

all: raytrace

raytrace: $(OBJ)/main.o $(OBJ)/parse.o $(OBJ)/raytrace.o $(OBJ)/image.o $(OBJ)/pixel.o $(OBJ)/EasyBMP.o
	$(CC) -std=c++11 $(DFLAGS) -o raytrace $(OBJ)/main.o $(OBJ)/parse.o $(OBJ)/raytrace.o $(OBJ)/image.o $(OBJ)/pixel.o $(OBJ)/EasyBMP.o

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