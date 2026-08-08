CXX ?= c++
CXXFLAGS ?= -O2 -std=c++14 -Wall -Wextra -pedantic
CPPFLAGS ?= -I/opt/homebrew/opt/eigen/include/eigen3

.PHONY: all clean hydrodynamics tensor stress

all: hydrodynamics piston tensor

hydrodynamics: particle_sim

particle_sim: main.cpp particle.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main.cpp -o particle_sim

piston: main_piston.cpp particle.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main_piston.cpp -o piston

tensor: matrix_base.o

matrix_base.o: matrix_base.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c matrix_base.cpp -o matrix_base.o

stress: stress_sim

stress_base.o: stress_base.cpp particle_base.cpp matrix_base.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c stress_base.cpp -o stress_base.o

stress_sim: stress_base.o
	$(CXX) $(CXXFLAGS) stress_base.o -o stress_sim

clean:
	rm -f particle_sim piston stress_sim *.o
