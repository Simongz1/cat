#include "ADComputeStressTensor.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>
#include "MathFVUtils.h"
#include <type_traits>
#include "MooseFunctorArguments.h"
#include "FaceInfo.h"

registerMooseObject("mlApp", ADComputeStressTensor);

InputParameters
ADComputeStressTensor::validParams(){

    InputParameters params = FunctorMaterial::validParams();
    params.addClassDescription("reconstructs the Lagrangian deformation gradient from the partial gradients of each component of the inverse map field");

    //request the full material point position vector
    // params.addRequiredParam<MooseFunctorName>("rhoX_x_name", "name of the material point position vector multiplied by the density in the x direction");
    // params.addRequiredParam<MooseFunctorName>("rhoX_y_name", "name of the material point position vector multiplied by the density in the y direction");
    // params.addRequiredParam<MooseFunctorName>("rhoX_z_name", "name of the material point position vector multiplied by the density in the z direction");
    params.addParam<MooseFunctorName>("X00", "X00", "name of the 00 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X01", "X01", "name of the 01 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X02", "X02", "name of the 02 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X10", "X10", "name of the 10 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X11", "X11", "name of the 11 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X12", "X12", "name of the 12 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X20", "X20", "name of the 20 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X21", "X21", "name of the 21 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X22", "X22", "name of the 22 component of the inverse deformation gradient");

    params.addRequiredParam<MooseFunctorName>("lambda_name", "name of the lambda lame constant");
    params.addRequiredParam<MooseFunctorName>("mu_name", "name of the mu lame constant");

    //also provide directly momentum components and density
    params.addRequiredParam<MooseFunctorName>("momentum_vector_name", "name of the momentum vector");
    params.addRequiredParam<MooseFunctorName>("density_name", "name of the density variable");

    ////
    params.addRequiredParam<MooseFunctorName>("stress_col1_name", "name of the stress column 1");
    params.addRequiredParam<MooseFunctorName>("stress_col2_name", "name of the stress column 2");
    params.addRequiredParam<MooseFunctorName>("stress_col3_name", "name of the stress column 3");

    params.addRequiredParam<MooseFunctorName>("T_row1_name", "name of the mechanical flux tensor component in the x direction");
    params.addRequiredParam<MooseFunctorName>("T_row2_name", "name of the mechanical flux tensor component in the y direction");
    params.addRequiredParam<MooseFunctorName>("T_row3_name", "name of the mechanical flux tensor component in the z direction");
    params.addParam<std::string>("strain_type", "small_strain", "type of strain formulation");

    return params;
}

ADComputeStressTensor::ADComputeStressTensor(const InputParameters &params)
    : FunctorMaterial(params),

    //obtain the 9 components
      _X00(getFunctor<ADReal>("X00")),
      _X01(getFunctor<ADReal>("X01")),
      _X02(getFunctor<ADReal>("X02")),
      _X10(getFunctor<ADReal>("X10")),
      _X11(getFunctor<ADReal>("X11")),
      _X12(getFunctor<ADReal>("X12")),
      _X20(getFunctor<ADReal>("X20")),
      _X21(getFunctor<ADReal>("X21")),
      _X22(getFunctor<ADReal>("X22")),
      _lambda(getFunctor<ADReal>("lambda_name")),
      _mu(getFunctor<ADReal>("mu_name")),
      _m(getFunctor<ADRealVectorValue>("momentum_vector_name")),
      _rho(getFunctor<ADReal>("density_name")),
      _strain_type(getParam<std::string>("strain_type"))
{   

    //for a lambda that creates stress internally
    auto computeFullStressTensor = [this](const auto & r, const auto & state) -> ADRankTwoTensor{
        //obtain individual components of the inverse deformation gradient
        const ADReal X00 = _X00(r, state);
        const ADReal X01 = _X01(r, state);
        const ADReal X02 = _X02(r, state);
        const ADReal X10 = _X10(r, state);
        const ADReal X11 = _X11(r, state);
        const ADReal X12 = _X12(r, state);
        const ADReal X20 = _X20(r, state);
        const ADReal X21 = _X21(r, state);
        const ADReal X22 = _X22(r, state);

        //assemble the full inverse deformation gradient
        ADRankTwoTensor X;
        X(0,0) = X00; X(0,1) = X01; X(0,2) = X02;
        X(1,0) = X10; X(1,1) = X11; X(1,2) = X12;
        X(2,0) = X20; X(2,1) = X21; X(2,2) = X22;

        //form actual deformation gradient
        ADRankTwoTensor I; I.setToIdentity();  

        //obtain constants
        const ADReal lambda = _lambda(r, state);
        const ADReal mu = _mu(r, state);

        ADRankTwoTensor stress;

        if (_strain_type == "small_strain"){
            ADRankTwoTensor epsilon = 0.5 * (X + X.transpose()) - I;

            //stress from small strain
            stress = lambda * epsilon.trace() * I + 2 * mu * epsilon;
        }
        else if (_strain_type == "large_strain"){

            //obtain deformation gradient
            ADRankTwoTensor F = X.inverse();

            //compute jacobian of the deformation gradient
            ADReal J = F.det();

            //compute right cauchy green deformation tensor
            ADRankTwoTensor Cinv = X * X.transpose();
            
            //compute second piola
            ADRankTwoTensor PK2;
            PK2 = lambda * MetaPhysicL::log(MetaPhysicL::max(J, 1e-6)) * Cinv + mu * (I - Cinv);

            //compute cauchy stress
            stress = (1. / J) * F * PK2 * F.transpose();
        }
        else {
            mooseError("Strain has to be either large_strain (defualt), or small strain !!");
        }
        return - stress;
    };

    addFunctorProperty<ADRealVectorValue>(getParam<MooseFunctorName>("stress_col1_name"),
        [this, computeFullStressTensor](const auto & r, const auto & state) -> ADRealVectorValue{
            //call the full stress tensor calculation
            const ADRankTwoTensor sigma = computeFullStressTensor(r, state);

            //define the vector
            ADRealVectorValue col1;
            for (unsigned int i = 0; i < 3; ++i){
                col1(i) = sigma(i, 0);
            }
            return col1;
        });

    addFunctorProperty<ADRealVectorValue>(getParam<MooseFunctorName>("stress_col2_name"),
        [this, computeFullStressTensor](const auto & r, const auto & state) -> ADRealVectorValue{
            //call the full stress tensor calculation
            const ADRankTwoTensor sigma = computeFullStressTensor(r, state);

            //define the vector
            ADRealVectorValue col2;
            for (unsigned int i = 0; i < 3; ++i){
                col2(i) = sigma(i, 1);
            }
            return col2;
        });

    addFunctorProperty<ADRealVectorValue>(getParam<MooseFunctorName>("stress_col3_name"),
        [this, computeFullStressTensor](const auto & r, const auto & state) -> ADRealVectorValue{
            //call the full stress tensor calculation
            const ADRankTwoTensor sigma = computeFullStressTensor(r, state);

            //define the vector
            ADRealVectorValue col3;
            for (unsigned int i = 0; i < 3; ++i){
                col3(i) = sigma(i, 2);
            }
            return col3;
        });

    //now compute the full mechanical flux tensor and then expose components
    
    //first, internally compute a full tensor
    auto computeMechanicalFluxTensor = [this, computeFullStressTensor](const auto & r, const auto & state) -> ADRankTwoTensor{
        //obtain the full stress tensor
        const ADRankTwoTensor sigma = computeFullStressTensor(r, state);

        //obtain the momentum vector
        const ADRealVectorValue m = _m(r, state);
        const ADReal rho = _rho(r, state);

        //compute the tensor

        ADRankTwoTensor T;
        for (unsigned int i = 0; i < 3; ++i){
            for (unsigned int j = 0; j < 3; ++j){
                T(i,j) = m(i) * m(j) / rho;
                T(i,j) -= sigma(i,j);
            }
        }
        return T;
    };

    //now declare the rows using this internal lambda
    addFunctorProperty<ADRealVectorValue>(getParam<MooseFunctorName>("T_row1_name"),
        [this, computeMechanicalFluxTensor](const auto & r, const auto & state) -> ADRealVectorValue{
            // //call the full mechanical flux tensor
            const ADRankTwoTensor T = computeMechanicalFluxTensor(r, state);

            //declare and populate first row
            ADRealVectorValue row1;
            for (unsigned int j = 0; j < 3; ++j){
                row1(j) = T(0,j);
            }
            return row1;
        });

    addFunctorProperty<ADRealVectorValue>(getParam<MooseFunctorName>("T_row2_name"),
        [this, computeMechanicalFluxTensor](const auto & r, const auto & state) -> ADRealVectorValue{
            //call the full mechanical flux tensor
            const ADRankTwoTensor T = computeMechanicalFluxTensor(r, state);

            //declare and populate first row
            ADRealVectorValue row2;
            for (unsigned int j = 0; j < 3; ++j){
                row2(j) = T(1,j);
            }
            return row2;
        });

    addFunctorProperty<ADRealVectorValue>(getParam<MooseFunctorName>("T_row3_name"),
        [this, computeMechanicalFluxTensor](const auto & r, const auto & state) -> ADRealVectorValue{
            //call the full mechanical flux tensor
            const ADRankTwoTensor T = computeMechanicalFluxTensor(r, state);

            //declare and populate first row
            ADRealVectorValue row3;
            for (unsigned int j = 0; j < 3; ++j){
                row3(j) = T(2,j);
            }
            return row3;
        });
}