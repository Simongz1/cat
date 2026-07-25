#include <iostream>
#include "particle.cpp"
#include <array>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

using Force = std::array<float, 2>;
constexpr float pi = 3.14159265358979323846f;

struct Lattice {
    std::vector<particle> particles;
    std::vector<int> grid_to_particle;
};

struct ModelParameters {
    float stiffness = 70.0f;
    float piston_velocity = -4.0f;
    float piston_ramp_time = 0.08f;
    bool use_roller_boundaries = true;
};

bool inside_center_hole(
    float x,
    float y,
    float xmin,
    float xmax,
    float ymin,
    float ymax,
    float hole_radius
){
    const float x_center = 0.5f * (xmin + xmax);
    const float y_center = 0.5f * (ymin + ymax);
    const float dx = x - x_center;
    const float dy = y - y_center;
    return std::sqrt(dx * dx + dy * dy) < hole_radius;
}

int original_i(const particle& p, int ny){
    return p.body_id / ny;
}

int original_j(const particle& p, int ny){
    return p.body_id % ny;
}

bool is_right_boundary_particle(const particle& p, int nx, int ny){
    return original_i(p, ny) == nx - 1;
}

bool is_left_boundary_particle(const particle& p, int ny){
    return original_i(p, ny) == 0;
}

bool is_top_or_bottom_particle(const particle& p, int ny){
    const int j = original_j(p, ny);
    return j == 0 || j == ny - 1;
}

Lattice initialize_particles(
    float xmin,
    float xmax,
    int nx,
    float ymin,
    float ymax,
    int ny,
    const std::string& geometry,
    float hole_radius
){
    Lattice lattice;
    std::vector<particle>& particles = lattice.particles;
    std::vector<int>& grid_to_particle = lattice.grid_to_particle;
    particles.reserve(nx * ny);
    grid_to_particle.assign(nx * ny, -1);

    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            particle p;
            p.x = nx > 1 ? xmin + (xmax - xmin) * i / (nx - 1) : xmin;
            p.y = ny > 1 ? ymin + (ymax - ymin) * j / (ny - 1) : ymin;

            if (geometry == "hole" && inside_center_hole(p.x, p.y, xmin, xmax, ymin, ymax, hole_radius)){
                continue;
            }

            p.x_old = p.x;
            p.y_old = p.y;
            p.x_ref = p.x;
            p.y_ref = p.y;
            p.body_id = i * ny + j;

            grid_to_particle[i * ny + j] = particles.size();
            particles.push_back(p);
        }
    }

    return lattice;
}

float axis_spacing(float min_value, float max_value, int npoints){
    return npoints > 1 ? (max_value - min_value) / (npoints - 1) : max_value - min_value;
}

float clamp_value(float value, float min_value, float max_value){
    if (value < min_value){
        return min_value;
    }
    if (value > max_value){
        return max_value;
    }
    return value;
}

void add_spring_force(
    std::vector<Force>& forces,
    const std::vector<particle>& particles,
    int a,
    int b,
    float rest_length,
    float stiffness
){
    constexpr float minimum_distance = 1.0e-3f;

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
    const std::vector<int>& grid_to_particle,
    int nx,
    int ny,
    float x_spacing,
    float y_spacing,
    const ModelParameters& model
){
    std::vector<Force> forces(particles.size(), Force{0.0f, 0.0f});

    if (particles.empty() || nx < 1 || ny < 1){
        return forces;
    }

    const float diagonal_spacing = std::sqrt(x_spacing * x_spacing + y_spacing * y_spacing);

    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            const int index = grid_to_particle[i * ny + j];
            if (index < 0){
                continue;
            }

            if (i + 1 < nx){
                const int neighbor = grid_to_particle[(i + 1) * ny + j];
                if (neighbor >= 0){
                    add_spring_force(forces, particles, index, neighbor, x_spacing, model.stiffness);
                }
            }
            if (j + 1 < ny){
                const int neighbor = grid_to_particle[i * ny + (j + 1)];
                if (neighbor >= 0){
                    add_spring_force(forces, particles, index, neighbor, y_spacing, model.stiffness);
                }
            }
            if (i + 1 < nx && j + 1 < ny){
                const int neighbor = grid_to_particle[(i + 1) * ny + (j + 1)];
                if (neighbor >= 0){
                    add_spring_force(forces, particles, index, neighbor, diagonal_spacing, model.stiffness);
                }
            }
            if (i + 1 < nx && j - 1 >= 0){
                const int neighbor = grid_to_particle[(i + 1) * ny + (j - 1)];
                if (neighbor >= 0){
                    add_spring_force(forces, particles, index, neighbor, diagonal_spacing, model.stiffness);
                }
            }
        }
    }

    return forces;
}

float ramped_piston_velocity(float time, const ModelParameters& model){
    if (model.piston_ramp_time <= 0.0f){
        return model.piston_velocity;
    }

    const float s = clamp_value(time / model.piston_ramp_time, 0.0f, 1.0f);
    return model.piston_velocity * (0.5f - 0.5f * std::cos(pi * s));
}

float piston_displacement(float time, const ModelParameters& model){
    if (time <= 0.0f){
        return 0.0f;
    }

    if (model.piston_ramp_time <= 0.0f){
        return model.piston_velocity * time;
    }

    if (time <= model.piston_ramp_time){
        const float omega = pi / model.piston_ramp_time;
        return model.piston_velocity * (0.5f * time - 0.5f * std::sin(omega * time) / omega);
    }

    const float ramp_displacement = piston_displacement(model.piston_ramp_time, model);
    return ramp_displacement + model.piston_velocity * (time - model.piston_ramp_time);
}

void apply_boundary_conditions(
    std::vector<particle>& particles,
    float time,
    int nx,
    int ny,
    const ModelParameters& model
){
    const float piston_velocity = ramped_piston_velocity(time, model);
    const float displacement = piston_displacement(time, model);

    for (particle& p : particles){
        if (is_right_boundary_particle(p, nx, ny)){
            p.x = p.x_ref + displacement;
            p.x_old = p.x;
            p.vx = piston_velocity;
            p.vx_old = p.vx;
            p.ax = 0.0f;
            p.ax_old = 0.0f;
        }

        if (is_left_boundary_particle(p, ny)){
            p.x = p.x_ref;
            p.x_old = p.x;
            p.vx = 0.0f;
            p.vx_old = 0.0f;
            p.ax = 0.0f;
            p.ax_old = 0.0f;
        }

        if (model.use_roller_boundaries && is_top_or_bottom_particle(p, ny)){
            p.y = p.y_ref;
            p.y_old = p.y;
            p.vy = 0.0f;
            p.vy_old = 0.0f;
            p.ay = 0.0f;
            p.ay_old = 0.0f;
        }
    }
}

void leapfrog(
    std::vector<particle>& particles,
    const std::vector<int>& grid_to_particle,
    float dt,
    float ymax,
    float ymin,
    float xmax,
    float xmin,
    int nx,
    int ny,
    float x_spacing,
    float y_spacing,
    float time,
    const ModelParameters& model,
    float& KE
){
    const std::vector<Force> forces = interactions(particles, grid_to_particle, nx, ny, x_spacing, y_spacing, model);

    for (unsigned int i = 0; i < particles.size(); ++i){
        particle& current_particle = particles[i];

        current_particle.ax = forces[i][0] / current_particle.m;
        current_particle.ay = forces[i][1] / current_particle.m;

        current_particle.vx = current_particle.vx_old + 0.5f * dt * (current_particle.ax + current_particle.ax_old);
        current_particle.vy = current_particle.vy_old + 0.5f * dt * (current_particle.ay + current_particle.ay_old);

        current_particle.x = current_particle.x_old + dt * current_particle.vx + 0.5f * current_particle.ax * dt * dt;
        current_particle.y = current_particle.y_old + dt * current_particle.vy + 0.5f * current_particle.ay * dt * dt;

        if (current_particle.y > ymax){
            current_particle.y = ymax - (current_particle.y - ymax);
            current_particle.vy *= -1.0f;
        }
        else if (current_particle.y < ymin){
            current_particle.y = ymin + (ymin - current_particle.y);
            current_particle.vy *= -1.0f;
        }

        if (current_particle.x > xmax){
            current_particle.x = xmax - (current_particle.x - xmax);
            current_particle.vx *= -1.0f;
        }
        else if (current_particle.x < xmin){
            current_particle.x = xmin + (xmin - current_particle.x);
            current_particle.vx *= -1.0f;
        }

        KE += 0.5f * current_particle.m * (current_particle.vx * current_particle.vx + current_particle.vy * current_particle.vy);

        current_particle.x_old = current_particle.x;
        current_particle.y_old = current_particle.y;
        current_particle.vx_old = current_particle.vx;
        current_particle.vy_old = current_particle.vy;
        current_particle.ax_old = current_particle.ax;
        current_particle.ay_old = current_particle.ay;
    }

    apply_boundary_conditions(particles, time + dt, nx, ny, model);

    const float avg_KE = KE / particles.size();
    std::cout << "Current kinetic energy " << avg_KE << " for " << particles.size() << " particles " << std::endl;
}

int main(int argc, char* argv[]){
    if (argc != 9 && argc != 10 && argc != 11){
        std::cerr << "Usage: " << argv[0]
                  << " xmin xmax nx ymin ymax ny dt nsteps [solid|hole] [hole_radius]" << std::endl;
        return 1;
    }

    const float xmin = std::stof(argv[1]);
    const float xmax = std::stof(argv[2]);
    const int nx = std::stoi(argv[3]);
    const float ymin = std::stof(argv[4]);
    const float ymax = std::stof(argv[5]);
    const int ny = std::stoi(argv[6]);
    const float dt = std::stof(argv[7]);
    const int nsteps = std::stoi(argv[8]);
    const std::string geometry = argc >= 10 ? argv[9] : "solid";

    const float x_spacing = axis_spacing(xmin, xmax, nx);
    const float y_spacing = axis_spacing(ymin, ymax, ny);
    const float default_hole_radius = 0.15f * std::min(xmax - xmin, ymax - ymin);
    const float hole_radius = argc == 11 ? std::stof(argv[10]) : default_hole_radius;

    if (geometry != "solid" && geometry != "hole"){
        std::cerr << "Unsupported geometry: " << geometry << std::endl;
        return 1;
    }

    const ModelParameters model;
    Lattice lattice = initialize_particles(xmin, xmax, nx, ymin, ymax, ny, geometry, hole_radius);
    std::vector<particle>& particles = lattice.particles;

    apply_boundary_conditions(particles, 0.0f, nx, ny, model);

    std::ofstream positions("positions_piston.csv");

    for (int t = 0; t < nsteps; ++t){
        std::cout << "STEP " << t << " EXECUTING INTEGRATION" << std::endl;
        float KE = 0.0f;
        const float time = t * dt;

        leapfrog(
            particles,
            lattice.grid_to_particle,
            dt,
            ymax,
            ymin,
            xmax,
            xmin,
            nx,
            ny,
            x_spacing,
            y_spacing,
            time,
            model,
            KE
        );

        if (t == 0){
            positions << "time,x,y,vmag,KE" << std::endl;
        }

        for (unsigned int i = 0; i < particles.size(); ++i){
            const float vx = particles[i].vx;
            const float vy = particles[i].vy;
            const float mag = std::sqrt(vx * vx + vy * vy);

            positions << time + dt << ","
                      << particles[i].x << ","
                      << particles[i].y << ","
                      << mag << ","
                      << KE / particles.size()
                      << std::endl;
        }
    }

    return 0;
}
