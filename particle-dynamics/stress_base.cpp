#include <iostream>
//call the particle definition
#include "particle_base.cpp"
#include <array>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

struct ModelParameters {
    bool use_prescribed_piston = false;
    float stiffness = 10.0f;
    float piston_velocity = -4.0f;
    float piston_ramp_time = 0.08f;
    bool use_roller_boundaries = true;
    float contact_stiffness = 100.0f;
    float contact_distance_factor = 1.0f;
};

//matrix based solver.
//there is an iterative computation of an interaction matrix
//that is converted into a linear system
//the interparticle force is computed using traction as \sigma \cdot n
//where n is the unit directional vector between each pair of particles
//within the radius of interaction

//define a function to compute interaction radius
void updateDistanceMatrix(
    const float interactionRadius,
    const std::vector<particle> &particles,
    Matrix<float, sizeof(particles), sizeof(particles)> &distanceMatrix,
    Matrix<float, sizeof(particles), sizeof(particles)> &distanceMatrixOld
){
    //obtain the current location for each particle
    for (unsigned int i = 0; i < sizeof(particles); ++i){
        for (unsigned int j = i + 1; j < sizeof(particles); ++j){
            //extract particle i and j
            particle particle_i = particles[i];
            particle particle_j = particles[j];

            //obtain the distance as the norm of the difference in current positions
            Matrix<float, 3 ,1> dx = particle_i.x - particle_j.x;

            //compute the norm
            float dx_norm = dx.norm();

            //save value in symmetric matrix
            distanceMatrix(i, j) = dx_norm;

            //set to zero particles that are outside of the interaction radius
            distanceMatrix(j, i) = distanceMatrix(i, j);

            //repeat for the old distance matrix

            Matrix<float, 3, 1> dx_old = particle_i.x_old - particle_j.x_old;
            float dx_old_norm = dx_old.norm();
            distanceMatrixOld(i, j) = dx_old_norm; 
            distanceMatrixOld(j, i) = distanceMatrixOld(i, j);

            //normalize the matrix by the interaction radius
            distanceMatrix(i, j) *= 1.0 / interactionRadius;
            distanceMatrixOld(i, j) *= 1.0 / interactionRadius;
        }
    }
}

//once we have updated the distance matrix, we now approximate the moment matrix

void updateWeights(
    const std::vector<particle> &particles,
    const Matrix<float, sizeof(particles), sizeof(particles)> &distanceMatrix,
    Matrix<float, sizeof(particles), sizeof(particles)> &weightsMatrix
){
    //compute weights
    for (unsigned int i = 0; i < sizeof(particles); ++i){
        for(unsigned int j = 0 ; j < sizeof(particles); ++j){
            weightsMatrix(i, j) = distanceMatrix(i, j) > 1.0f ? std::pow(1.0 - distanceMatrix(i, j), 4) * (1.0 + 4.0 * distanceMatrix(i, j)) : 0.0f;
        }
    }
}

//now we can define a function that computes the momentum and cross config tensors

void updateMomentumTensor(
    const std::vector<particle> &particles,
    const Matrix<float, sizeof(particles), sizeof(particles)> &momentumTensor,
    const Matrix<float, sizeof(particles), sizeof(particles)> &weightsMatrix,
    const Matrix<float, sizeof(particles), sizeof(particles)> &distanceMatrixOld
){
    //compute the momentum tensor for each particle
    for (unsigned int i = 0; i < sizeof(particles); ++i){
        particle particle_i = particles[i];
        Matrix<float, 3, 3> M_i = particle_i.M;
        for (unsigned int j = i + 1; j < sizeof(particles); ++j){
            //extract particles
            
            particle particle_j = particles[j];

            //update the momentum tensor
            
            M_i += weightsMatrix(i, j) * distanceMatrixOld(i, j) * distanceMatrixOld(j, i);
        }
    }
}

int main(){
    return 0;
}