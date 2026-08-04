// define a particle structure class
#include <vector>
#include "matrix_base.cpp"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct particle {
    /////////
    float x1 = 0.0, x2 = 0.0, x3 = 0.0;
    float x1_old = 0.0, x2_old = 0.0, x3_old = 0.0;
    Matrix<float, 3, 1> x = {x1, x2, x3};
    Matrix<float, 3, 1> x_old = {x1_old, x2_old, x3_old};
    /////////
    float v1 = 0.0, v2 = 0.0, v3 = 0.0;
    float v1_old = 0.0, v2_old = 0.0, v3_old = 0.0;
    Matrix<float, 3, 1> v = {v1, v2, v3};
    Matrix<float, 3, 1> v_old = {v1_old, v2_old, v3_old};
    /////////
    float a1 = 0.0, a2 = 0.0, a3 = 0.0;
    float a1_old = 0.0, a2_old = 0.0, a3_old = 0.0;
    Matrix<float, 3, 1> a = {a1, a2, a3};
    Matrix<float, 3, 1> a_old = {a1_old, a2_old, a3_old};
    /////////
    Matrix<float, 3, 3> F;
    float J = 1.0f;
    float m = 1000e-6;

    //initialize momentum tensor
    Matrix<float, 3, 3> M;
    Matrix<float, 3, 3> B;
};