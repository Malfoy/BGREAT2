# CC=/usr/bin/g++
CC=g++
#~ CC=clang++
CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors -fopenmp
LDFLAGS=-flto -lpthread -fopenmp


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4   -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=bgreat

all: $(EXEC)



aligner.o: aligner.cpp aligner.h utils.h alignerGreedy.cpp alignerExhaustive.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: utils.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

bgreat: bgreat.o   aligner.o utils.o
	$(CC) -o $@ $^ $(LDFLAGS)

getBigUnitigs.o: getBigUnitigs.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

bgreat.o: bgreat.cpp  aligner.h
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
