CC=g++
OS := $(shell uname)
ifeq ($(OS),Darwin)
  # Run MacOS commands
  CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors -lz
  LDFLAGS=-flto -lpthread -lz
else
  # check for Linux and add use of OpenMP
  CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors -fopenmp -lz
  LDFLAGS=-flto -lpthread -fopenmp -lz
endif


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4   -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=bgreat numbersToSequences sortPaths

all: $(EXEC)

sortPaths: sortPaths.o
	$(CC) -o $@ $^ $(LDFLAGS)

sortPaths.o: sortPaths.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

aligner.o: aligner.cpp aligner.h utils.h alignerGreedy.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: utils.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

bgreat: bgreat.o   aligner.o utils.o
	$(CC) -o $@ $^ $(LDFLAGS)

bgreat.o: bgreat.cpp  aligner.h
	$(CC) -o $@ -c $< $(CFLAGS)

numbersToSequences: numbersToSequences.o
	$(CC) -o $@ $^ $(LDFLAGS)

numbersToSequences.o: numbersToSequences.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
