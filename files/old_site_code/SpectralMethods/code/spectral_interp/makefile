CXX = g++
CXXFLAGS = -std=c++11 -O3 -march=native
LIBS = -lfftw3 -lm
#LIBFLAGS = 
DEBUG = -Wall -Werror

all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))
%.out: %.cpp makefile
	$(CXX) $(CXXFLAGS) $(LIBS) $(LIBFLAGS) $(DEBUG) $< -o $@ 
clean:
	rm -f *.out

.PHONY: all, clean
