// define a particle structure class
#include <vector>
#include "matrix_base.cpp"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct particle {
    /////////
    Matrix<float, 3, 1> x = {0, 0, 0};
    Matrix<float, 3, 1> x_old = {0, 0, 0};
    float x1 = x[0], x2 = x[1], x3 = x[2];
    float x1_old = x_old[0], x2_old = x_old[1], x3_old = x_old[2];
    /////////
    Matrix<float, 3, 1> v = {0, 0, 0};
    Matrix<float, 3, 1> v_old = {0, 0, 0};
    float v1 = v[0], v2 = v[1], v3 = v[2];
    float v1_old = v_old[0], v2_old = v_old[1], v3_old = v_old [2];
    /////////
    Matrix<float, 3, 1> a = {0, 0, 0};
    Matrix<float, 3, 1> a_old = {0, 0, 0};
    float a1 = a[0], a2 = a[1], a3 = a[2];
    float a1_old = a_old[0], a2_old = a_old[1], a3_old = a_old[2];
    /////////
    Matrix<float, 3, 3> F;
    Matrix<float, 3, 3> F_old;
    Matrix<float, 3, 3> f;
    float J = 1.0f;
    float m = 1000e-6;

    //initialize momentum tensor
    Matrix<float, 3, 3> M;
    Matrix<float, 3, 3> B;

    static void initKinematics(Matrix<float, 3, 3> &F, Matrix<float, 3, 3> &F_old){
        for (int i = 0; i < 3; ++i){
            F(i, i) = 1;
            F_old(i, i) = 1;
        }
    };

};