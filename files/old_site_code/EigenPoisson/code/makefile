CXX = g++
CXXFLAGS = -std=c++11 -O3 -march=native
LIBS = -I /usr/include/eigen3/
#LIBFLAGS = -DEIGEN_DONT_PARALLELIZE 
DEBUG = -Wall -Werror

all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))
%.out: %.cpp makefile
	$(CXX) $(CXXFLAGS) $(LIBS) $(LIBFLAGS) $(DEBUG) $< -o $@ -lm
clean:
	rm -f *.out

.PHONY: all, clean
