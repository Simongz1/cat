#include <iostream>
//call the particle definition
#include "particle.cpp"
#include <array>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unordered_map>
using namespace std;

constexpr float pi = 3.14159265358979323846f;

//define a particle initialization function
//this function takes as input the grid parameters and outputs a container of particles

struct Lattice {
    std::vector<particle> particles;
    std::vector<int> grid_to_particle;
};

struct ModelParameters {
    bool use_prescribed_piston = false;
    float stiffness = 10.0f;
    float piston_velocity = -4.0f;
    float piston_ramp_time = 0.08f;
    bool use_roller_boundaries = true;
    float contact_stiffness = 100.0f;
    float contact_distance_factor = 1.0f;
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
    const float x_center = 0.7f * (xmin + xmax);
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

            if (i == nx - 1){
                p.vx = -1000.0f;
                p.vx_old = p.vx;
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
    constexpr float minimum_distance = 1e-6;

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

struct Cell {
    int x;
    int y;

    bool operator==(const Cell& other) const {
        return x == other.x && y == other.y;
    }
};

struct CellHash {
    std::size_t operator()(const Cell& cell) const {
        const std::size_t hx = std::hash<int>{}(cell.x);
        const std::size_t hy = std::hash<int>{}(cell.y);
        return hx ^ (hy + 0x9e3779b9U + (hx << 6) + (hx >> 2));
    }
};

void add_contact_force(
    std::vector<Force>& forces,
    const std::vector<particle>& particles,
    int a,
    int b,
    float contact_distance,
    float contact_stiffness
){
    float dx = particles[a].x - particles[b].x;
    float dy = particles[a].y - particles[b].y;
    const float distance = std::sqrt(dx * dx + dy * dy);

    if (distance >= contact_distance){
        return;
    }

    // At exact overlap, use the particles' original separation to define
    // a stable direction instead of silently dropping the contact force.
    constexpr float minimum_distance = 1e-8f;
    float direction_length = distance;
    if (direction_length <= minimum_distance){
        dx = particles[a].x_ref - particles[b].x_ref;
        dy = particles[a].y_ref - particles[b].y_ref;
        direction_length = std::sqrt(dx * dx + dy * dy);
        if (direction_length <= minimum_distance){
            dx = 1.0f;
            dy = 0.0f;
            direction_length = 1.0f;
        }
    }

    const float magnitude =
        contact_stiffness * (contact_distance - distance);
    const float fx = magnitude * dx / direction_length;
    const float fy = magnitude * dy / direction_length;

    forces[a][0] += fx;
    forces[a][1] += fy;
    forces[b][0] -= fx;
    forces[b][1] -= fy;
}

void add_particle_contacts(
    std::vector<Force>& forces,
    const std::vector<particle>& particles,
    float contact_distance,
    float contact_stiffness
){
    if (contact_distance <= 0.0f || contact_stiffness <= 0.0f){
        return;
    }

    // Each cell is one cutoff wide, so possible contacts can only occur in
    // the particle's own cell or one of its eight neighboring cells.
    std::unordered_map<Cell, std::vector<int>, CellHash> cells;
    cells.reserve(particles.size());

    for (int index = 0; index < static_cast<int>(particles.size()); ++index){
        const Cell cell{
            static_cast<int>(std::floor(particles[index].x / contact_distance)),
            static_cast<int>(std::floor(particles[index].y / contact_distance))
        };

        for (int offset_x = -1; offset_x <= 1; ++offset_x){
            for (int offset_y = -1; offset_y <= 1; ++offset_y){
                const Cell neighbor_cell{
                    cell.x + offset_x,
                    cell.y + offset_y
                };
                const auto neighbor = cells.find(neighbor_cell);
                if (neighbor == cells.end()){
                    continue;
                }

                for (const int other : neighbor->second){
                    add_contact_force(
                        forces,
                        particles,
                        index,
                        other,
                        contact_distance,
                        contact_stiffness
                    );
                }
            }
        }

        cells[cell].push_back(index);
    }
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

    const float contact_distance =
        model.contact_distance_factor * std::min(x_spacing, y_spacing);
    add_particle_contacts(
        forces,
        particles,
        contact_distance,
        model.contact_stiffness
    );

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
    if (!model.use_prescribed_piston){
        return;
    }

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

//define a simple leap frog integration for particle motion

void leapfrog(
    std::vector<particle> &particles,
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
    float &KE
){
    //evaliate acceleration
    const std::vector<Force> forces = interactions(particles, grid_to_particle, nx, ny, x_spacing, y_spacing, model);

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
        if (current_particle.y >= ymax){
            current_particle.y = ymax - (current_particle.y - ymax);
            current_particle.vy *= -1.0f;
        }
        else if (current_particle.y <= ymin){
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

    apply_boundary_conditions(particles, time + dt, nx, ny, model);

    float avg_KE = KE / particles.size();
    std::cout << "Current kinetic energy " << avg_KE << " for " << particles.size() << " particles " << std::endl;

}

int main(int argc, char* argv[]){
    if (argc < 9 || argc > 14){
        std::cerr << "Usage: " << argv[0]
                  << " xmin xmax nx ymin ymax ny dt nsteps [solid|hole]"
                  << " [hole_radius] [free|piston] [piston_velocity] [piston_ramp_time]" << std::endl;
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
    int next_arg = 9;
    std::string geometry = "solid";
    if (argc > next_arg){
        geometry = argv[next_arg];
        ++next_arg;
    }

    std::cout << "Domain x-range: " << xmin << ", " << xmax << std::endl;
    std::cout << "Domain y-range: " << ymin << ", " << ymax << std::endl;
    std::cout << "Running " << nsteps << " integration steps with dt = " << dt << std::endl;

    const float x_spacing = axis_spacing(xmin, xmax, nx);
    const float y_spacing = axis_spacing(ymin, ymax, ny);
    const float default_hole_radius = 0.15f * std::min(xmax - xmin, ymax - ymin);
    float hole_radius = default_hole_radius;
    if (geometry == "hole" && argc > next_arg){
        const std::string candidate = argv[next_arg];
        if (candidate != "free" && candidate != "piston"){
            hole_radius = std::stof(candidate);
            ++next_arg;
        }
    }

    ModelParameters model;
    std::string boundary_mode = "free";
    if (argc > next_arg){
        boundary_mode = argv[next_arg];
        ++next_arg;
    }

    if (boundary_mode == "piston"){
        model.use_prescribed_piston = true;
        if (argc > next_arg){
            model.piston_velocity = std::stof(argv[next_arg]);
            ++next_arg;
        }
        if (argc > next_arg){
            model.piston_ramp_time = std::stof(argv[next_arg]);
            ++next_arg;
        }
    }

    std::cout << "Particle spacing x: " << x_spacing << std::endl;
    std::cout << "Particle spacing y: " << y_spacing << std::endl;
    std::cout << "Geometry: " << geometry << std::endl;
    if (geometry == "hole"){
        std::cout << "Hole radius: " << hole_radius << std::endl;
    }

    if (geometry != "solid" && geometry != "hole"){
        std::cerr << "Unsupported geometry: " << geometry << std::endl;
        return 1;
    }

    if (boundary_mode != "free" && boundary_mode != "piston"){
        std::cerr << "Unsupported boundary mode: " << boundary_mode << std::endl;
        return 1;
    }

    if (next_arg != argc){
        std::cerr << "Too many arguments provided." << std::endl;
        return 1;
    }

    std::cout << "Boundary mode: " << boundary_mode << std::endl;
    if (model.use_prescribed_piston){
        std::cout << "Piston velocity: " << model.piston_velocity << std::endl;
        std::cout << "Piston ramp time: " << model.piston_ramp_time << std::endl;
    }

    Lattice lattice = initialize_particles(xmin, xmax, nx, ymin, ymax, ny, geometry, hole_radius);
    std::vector<particle>& particles = lattice.particles;

    std::cout << "Initialized " << particles.size() << " particles" << std::endl;

    if (model.use_prescribed_piston){
        for (particle& p : particles){
            if (is_right_boundary_particle(p, nx, ny)){
                p.vx = 0.0f;
                p.vx_old = 0.0f;
            }
        }
        apply_boundary_conditions(particles, 0.0f, nx, ny, model);
    }

    //generate file container
    std::ofstream positions;
    positions.open("positions.csv");

    //iterate over time steps
    for (int t = 0; t < nsteps; ++t){
        //plot information to console

        std::cout << "STEP " << t << " EXECUTING INTEGRATION" << std::endl;
        float KE = 0.0;
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

            positions << time + dt << "," << particles[i].x << "," << particles[i].y << "," << mag << "," << KE / particles.size() << std::endl;
        }
    }
    positions.close();

    return 0;
}
