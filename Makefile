CXX ?= c++
CXXFLAGS ?= -O2 -std=c++14 -Wall -Wextra -pedantic

.PHONY: all clean

all: hydrodynamics piston

hydrodynamics: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o hydrodynamics

piston: main_piston.cpp
	$(CXX) $(CXXFLAGS) main_piston.cpp -o piston

clean:
	rm -f hydrodynamics piston *.o
