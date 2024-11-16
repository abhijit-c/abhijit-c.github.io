CXX = gcc
CXXFLAGS = -O3 -march=native
LIBS = -fopenmp 
#LIBFLAGS = -DEIGEN_DONT_PARALLELIZE 
DEBUG = -Wall 

all: $(patsubst %.c, %.out, $(wildcard *.c))
%.out: %.c makefile
	$(CXX) $(CXXFLAGS) $(DEBUG) $< -o $@ -lm -fopenmp
clean:
	rm -f *.out
	rm -f *.pgm
	rm -f *.png

.PHONY: all, clean
