CXX ?= c++
CXXFLAGS ?= -O2 -std=c++14 -Wall -Wextra -pedantic

.PHONY: all clean hydrodynamics

all: hydrodynamics piston

hydrodynamics: particle_sim

particle_sim: main.cpp particle.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o particle_sim

piston: main_piston.cpp particle.cpp
	$(CXX) $(CXXFLAGS) main_piston.cpp -o piston

clean:
	rm -f particle_sim piston *.o
