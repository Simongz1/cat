#include <iostream>
#include <matplotlibcpp.h>
//call the particle definition
#include "particle.cpp"
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

int main(){

    //define domain size
    float xmin = 0.;
    float xmax = 10.;
    int nx = 10;

    float ymin = 0.;
    float ymax = 10.;
    int ny = 10;
    
    //define a particle container
    std::vector<particle> particles;

    //define particles
    for (unsigned int i = 0; i < 10; ++i){
        for (unsigned int j = 0; j < 10; ++j){
            particle p;
            p.x = xmin + (xmax - xmin) * i / nx;
            p.y = ymin + (ymax - ymin) * j / ny;

            particles.push_back(p);
        }
    }

    matplot::scatter({particles[0].x}, {particles[0].y});
    matplot::show();
    return 0;
}