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
#include <Eigen/Sparse>
#include "geometry.cpp"
using namespace std;
using namespace Eigen;
using SpIndex = Eigen::Index;

//sparse matrix alias

using ParticleSparseMatrix = Eigen::SparseMatrix<float, Eigen::RowMajor>;

struct ModelParameters {
    bool use_prescribed_piston = false;
    float stiffness = 10.0f;
    float piston_velocity = -4.0f;
    float piston_ramp_time = 0.08f;
    bool use_roller_boundaries = true;
    float contact_stiffness = 100.0f;
    float contact_distance_factor = 1.0f;
};

struct SimulationParameters {
    const float dt = 1e-2;
    const float ti = 0.0;
    const float tf = 100.0;
    const float interactionRadius = 0.05;
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
    std::vector<particle> &particles,
    ParticleSparseMatrix &distanceMatrix,
    ParticleSparseMatrix &distanceMatrixOld
){  
    for (unsigned int i = 0; i < particles.size(); ++i){
        for (unsigned int j = i + 1; j < particles.size(); ++j){
            //define current and old entries
            std::vector<Eigen::Triplet<float>> currentInteracting;
            std::vector<Eigen::Triplet<float>> oldInteracting;

            //compute the actual differences
            Eigen::Vector3f currentDifference = particles[j].x - particles[i].x;
            Eigen::Vector3f oldDifference = particles[j].x_old - particles[i].x_old;

            //compute norms
            float currentDistance = currentDifference.norm();
            float oldDistance = oldDifference.norm();

            //check for interaction
            if (currentDistance < interactionRadius){
                //normalize
                float normalizedCurrentDistance = currentDistance / interactionRadius;
                float normalizedOldDistance = oldDistance / interactionRadius;

                //append to interacting list
                currentInteracting.emplace_back(i, j, normalizedCurrentDistance);
                currentInteracting.emplace_back(j, i, normalizedCurrentDistance);

                oldInteracting.emplace_back(i, j, normalizedOldDistance);
                oldInteracting.emplace_back(j, i, normalizedOldDistance);

                std::cout << "Particles " << i << " and " << j << " are interacting";
            }

            //set values on sparse matrix
            distanceMatrix.setFromTriplets(currentInteracting.begin(), currentInteracting.end());
            distanceMatrixOld.setFromTriplets(oldInteracting.begin(), oldInteracting.end());
        }
    }
}

//once we have updated the distance matrix, we now approximate the moment matrix

void updateWeights(
    std::vector<particle> &particles,
    const ParticleSparseMatrix &distanceMatrix,
    ParticleSparseMatrix &weightsMatrix,
    ParticleSparseMatrix &distanceMatrixOld
){
    //initialize sparse weight entries
    std::vector<Eigen::Triplet<float>> weightEntries;

    //iterate over sparse interactions
    for (Eigen::Index i = 0; i < distanceMatrixOld.outerSize(); ++i){
        for (ParticleSparseMatrix::InnerIterator entry(distanceMatrixOld, i); entry; ++entry){
            //form j index
            Eigen::Index j = entry.col();
            float q = entry.value();

            //evaluate distance
            if (q < 1.0f){
                //form distance kernel
                float m = 1.0 - q;
                float w = std::pow(m, 4.0) * (1.0 + 4.0 * m);
                
                //append to the tripled the location and weight values
                weightEntries.emplace_back(i, j, w);
            }
        }
    }

    weightsMatrix.setFromTriplets(weightEntries.begin(), weightEntries.end());
}

//now we can define a function that computes the momentum and cross config tensors

void updateMomentumTensors(
    std::vector<particle> &particles,
    const ParticleSparseMatrix &weightsMatrix,
    const ParticleSparseMatrix &distanceMatrixOld,
    const ParticleSparseMatrix &distanceMatrix
){
    //compute the momentum tensor for each particle
    //for sparse matrices, we only iterate through recorded index values

    //sparse indexes
    for (SpIndex i = 0; i < weightsMatrix.outerSize(); ++i){

        //obtain particle here
        particle currentParticle = particles[i];

        for (ParticleSparseMatrix::InnerIterator entry(weightsMatrix, i); entry; ++entry){

            //retrieve i, j weight
            SpIndex j = entry.col();
            float currentWeight = entry.value();

            //now that we have a nonzero pair, we compute the M and B tensors for particle i and j
            for (unsigned int k = 0; k < 3; ++k){
                for (unsigned int l = 0; l < 3; ++l){
                    currentParticle.M(k, l) += currentWeight * distanceMatrixOld.coeff(i, j) * distanceMatrixOld.coeff(j, i);
                    currentParticle.B(k, l) += currentWeight * distanceMatrix.coeff(i, j) * distanceMatrixOld.coeff(j, i);
                }
            }
        }
    }
}

//compute the incremental deformation gradient
void updateIncrementalDeformationGradient(
    std::vector<particle> &particles
){
    for (unsigned int i = 0; i < particles.size(); ++i){
        particles[i].f = particles[i].B * particles[i].M.inverse();
        particles[i].F = particles[i].f * particles[i].F_old;
    }
}

//once we have the deformation gradient we can compute the 

int main(){
    //we initially set initial structures
    GeometricParameters geom;
    std::vector<particle> particles = initParticles(geom);

    //std::cout << 60 * "=" << std::endl;
    std::cout << "Generated " << particles.size() << " particles" << std::endl;
    //std::cout << 60 * "=" << std::endl;

    //create an eigen index
    Eigen::Index particleCount = static_cast<Eigen::Index>(particles.size());

    //create the matrices
    ParticleSparseMatrix distanceMatrix(particleCount, particleCount);
    ParticleSparseMatrix distanceMatrixOld(particleCount, particleCount);
    ParticleSparseMatrix weightsMatrix(particleCount, particleCount);

    for (unsigned int step = 0; step < int(SimulationParameters().tf / SimulationParameters().dt); ++step){
        updateDistanceMatrix(0.01, particles, distanceMatrix, distanceMatrixOld);
        std::cout << "STEP " << step << std::endl;
    }

    return 0;
}