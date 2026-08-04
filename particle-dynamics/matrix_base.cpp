//this object declares a base class to perform
//templated matrix operations

#include <vector>
#include <cstddef>
#pragma once

// Three-dimensional vector.  Do not call this type `vector`, because that
// collides with std::vector in translation units that use this definition.
//define a class for tensors

//define a vector class with operations

class tensor {
    private:
        //define dimension and base tensor
        std::size_t dim;
        std::vector<std::vector<float>> A;

        static constexpr int item(int i, int j, int k){
            return ((i - j) * (j - k) * (k - i)) / 2;
        }

    public:
        //define 
        explicit tensor(std::size_t dimension)
            : dim(dimension),
              A(dimension, std::vector<float>(dimension, 0.0f))
        
        {}

        static tensor identity(std::size_t dimension){
            tensor result(dimension);

            for (std::size_t i = 0; i < dimension; ++i){
                result.A[i][i] = 1.0f;
            }

            return result;
        }

        std::vector<std::vector<float>> & values(){
            return A;
        }

        std::size_t dimension() const {
            return dim;
        }

        float det(std::size_t dimension) const {
            float det = 0.0f;
            //form the levi-civita symbol
            std::vector<std::vector<std::vector<float>>> epsilon(dimension);

            //define indexed expression
            if (dimension == 3){
                //levi civita tensor of rank 3
                for (unsigned int i = 0; i < dimension; ++i){
                    for (unsigned int j = 0; j < dimension; ++j){
                        for (unsigned int k = 0; k < dimension; ++k){
                            //define cases
                            epsilon[i][j][k] = item(i, j, k);
                        }
                    }
                }

                for (unsigned int i = 0; i < dimension; ++i){
                    for (unsigned int j = 0; j < dimension; ++j){
                        for (unsigned int k = 0; k < dimension; ++k){
                            det = epsilon[i][j][k] * A[1][i] * A[2][j] * A[3][k];
                        }
                    }
                }
            }
            return det;
        }

        float norm(std::vector<float> &v1, std::size_t dimension){
            float norm = 0.0f;

            for (unsigned int i = 0; i < dimension; ++i){
                norm += v1[i] * v1[i];
            }
            return std::sqrt(norm);
        }
};