#include <iostream>
//call the particle definition
#include "particle.cpp"
#include <array>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
using namespace std;

//define a particle initialization function
//this function takes as input the grid parameters and outputs a container of particles

std::vector<particle> initialize_particles(float xmin, float xmax, int nx, float ymin, float ymax, int ny){
    std::vector<particle> particles;
    particles.reserve(nx * ny);

    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            particle p;
            p.x = nx > 1 ? xmin + (xmax - xmin) * i / (nx - 1) : xmin;
            p.y = ny > 1 ? ymin + (ymax - ymin) * j / (ny - 1) : ymin;

            if (i == nx - 1){
                p.vx = -1.0f;
                p.vx_old = p.vx;
            }
            p.x_old = p.x;
            p.y_old = p.y;

            particles.push_back(p);
        }
    }

    return particles;
}

float axis_spacing(float min_value, float max_value, int npoints){
    return npoints > 1 ? (max_value - min_value) / (npoints - 1) : max_value - min_value;
}

//define a simple repulsion function for interparticle forces

using Force = std::array<float, 2>;

void add_spring_force(
    std::vector<Force>& forces,
    const std::vector<particle>& particles,
    int a,
    int b,
    float rest_length,
    float stiffness
){
    constexpr float minimum_distance = 1.0e-6f;

    const float dx = particles[a].x - particles[b].x;
    const float dy = particles[a].y - particles[b].y;
    const float distance = std::sqrt(dx * dx + dy * dy);

    if (distance <= minimum_distance){
        return;
    }

    const float extension = distance - rest_length;
    const float magnitude = -stiffness * extension;
    const float fx = magnitude * dx / distance;
    const float fy = magnitude * dy / distance;

    forces[a][0] += fx;
    forces[a][1] += fy;
    forces[b][0] -= fx;
    forces[b][1] -= fy;
}

std::vector<Force> interactions(
    const std::vector<particle>& particles,
    int nx,
    int ny,
    float x_spacing,
    float y_spacing
){
    constexpr float stiffness = 50.0f;
    std::vector<Force> forces(particles.size(), Force{0.0f, 0.0f});

    if (particles.empty() || nx < 1 || ny < 1){
        return forces;
    }

    const float diagonal_spacing = std::sqrt(x_spacing * x_spacing + y_spacing * y_spacing);

    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            const int index = i * ny + j;

            if (i + 1 < nx){
                add_spring_force(forces, particles, index, (i + 1) * ny + j, x_spacing, stiffness);
            }
            if (j + 1 < ny){
                add_spring_force(forces, particles, index, i * ny + (j + 1), y_spacing, stiffness);
            }
            if (i + 1 < nx && j + 1 < ny){
                add_spring_force(forces, particles, index, (i + 1) * ny + (j + 1), diagonal_spacing, stiffness);
            }
            if (i + 1 < nx && j - 1 >= 0){
                add_spring_force(forces, particles, index, (i + 1) * ny + (j - 1), diagonal_spacing, stiffness);
            }
        }
    }

    return forces;
}

//define a simple leap frog integration for particle motion

void leapfrog(
    std::vector<particle> &particles,
    float dt,
    float ymax,
    float ymin,
    float xmax,
    float xmin,
    int nx,
    int ny,
    float x_spacing,
    float y_spacing,
    float &KE
){
    //evaliate acceleration
    const std::vector<Force> forces = interactions(particles, nx, ny, x_spacing, y_spacing);

    //compute acceleration as F / m = a
    for (unsigned int i = 0; i < particles.size(); ++i){
        particle& current_particle = particles[i];

        current_particle.ax = forces[i][0] / current_particle.m;
        current_particle.ay = forces[i][1] / current_particle.m;

        //update velocity
        current_particle.vx = current_particle.vx_old + 0.5f * dt * (current_particle.ax + current_particle.ax_old);
        current_particle.vy = current_particle.vy_old + 0.5f * dt * (current_particle.ay + current_particle.ay_old);

        //update position
        current_particle.x = current_particle.x_old + dt * current_particle.vx + 0.5f * current_particle.ax * dt * dt;
        current_particle.y = current_particle.y_old + dt * current_particle.vy + 0.5f * current_particle.ay * dt * dt;

        // Reflect particles that cross the top or bottom boundary.
        if (current_particle.y > ymax){
            current_particle.y = ymax - (current_particle.y - ymax);
            current_particle.vy *= -1.0f;
        }
        else if (current_particle.y < ymin){
            current_particle.y = ymin + (ymin - current_particle.y);
            current_particle.vy *= -1.0f;
        }

        // Reflect particles that cross the left and right boundary.
        if (current_particle.x > xmax){
            current_particle.x = xmax - (current_particle.x - xmax);
            current_particle.vx *= -1.0f;
        }
        else if (current_particle.x < xmin){
            current_particle.x = xmin + (xmin - current_particle.x);
            current_particle.vx *= -1.0f;
        }

        //compute the total kinetic energy

        KE += 0.5f * current_particle.m * (current_particle.vx * current_particle.vx + current_particle.vy * current_particle.vy);
        
        //update old values
        current_particle.x_old = current_particle.x;
        current_particle.y_old = current_particle.y;

        current_particle.vx_old = current_particle.vx;
        current_particle.vy_old = current_particle.vy;

        current_particle.ax_old = current_particle.ax;
        current_particle.ay_old = current_particle.ay;
    }

    float avg_KE = KE / particles.size();
    std::cout << "Current kinetic energy " << avg_KE << " for " << particles.size() << " particles " << std::endl;

}

int main(int argc, char* argv[]){
    if (argc != 9){
        std::cerr << "Usage: " << argv[0]
                  << " xmin xmax nx ymin ymax ny dt nsteps" << std::endl;
        return 1;
    }

    //define domain size
    const float xmin = std::stof(argv[1]);
    const float xmax = std::stof(argv[2]);
    const int nx = std::stoi(argv[3]);

    const float ymin = std::stof(argv[4]);
    const float ymax = std::stof(argv[5]);
    const int ny = std::stoi(argv[6]);

    const float dt = std::stof(argv[7]);
    const int nsteps = std::stoi(argv[8]);

    std::cout << "Domain x-range: " << xmin << ", " << xmax << std::endl;
    std::cout << "Domain y-range: " << ymin << ", " << ymax << std::endl;
    std::cout << "Running " << nsteps << " integration steps with dt = " << dt << std::endl;

    const float x_spacing = axis_spacing(xmin, xmax, nx);
    const float y_spacing = axis_spacing(ymin, ymax, ny);
    std::cout << "Particle spacing x: " << x_spacing << std::endl;
    std::cout << "Particle spacing y: " << y_spacing << std::endl;

    std::vector<particle> particles = initialize_particles(xmin, xmax, nx, ymin, ymax, ny);

    //generate file container
    std::ofstream positions;
    positions.open("positions.csv");

    //iterate over time steps
    for (int t = 0; t < nsteps; ++t){
        //plot information to console

        std::cout << "STEP " << t << " EXECUTING INTEGRATION" << std::endl;
        float KE = 0.0;

        leapfrog(particles, dt, ymax, ymin, xmax, xmin, nx, ny, x_spacing, y_spacing, KE);
        //append column names as first row
        if (t == 0){
            positions << "time,x,y,vmag,KE" << std::endl;
        }

        //save locations on file
        for (unsigned int i = 0; i < particles.size(); ++i){
            //compute velocity magnitude
            float vx = particles[i].vx;
            float vy = particles[i].vy;
            float mag = std::sqrt(vx * vx + vy * vy);

            positions << t * dt << "," << particles[i].x << "," << particles[i].y << "," << mag << "," << KE / particles.size() << std::endl;
        }
    }
    positions.close();

    return 0;
}
