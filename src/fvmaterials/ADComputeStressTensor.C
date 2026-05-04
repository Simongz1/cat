#include "ADComputeStressTensor.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeStressTensor);

InputParameters
ADComputeStressTensor::validParams(){

    InputParameters params = FunctorMaterial::validParams();
    params.addClassDescription("reconstructs the Lagrangian deformation gradient from the partial gradients of each component of the inverse map field");
    params.addRequiredParam<MooseFunctorName>("map_x_gradient", "name of gradient of the x component of the inverse map");
    params.addRequiredParam<MooseFunctorName>("map_y_gradient", "name of gradient of the x component of the inverse map");
    params.addRequiredParam<MooseFunctorName>("map_z_gradient", "name of gradient of the x component of the inverse map");
    params.addRequiredParam<MooseFunctorName>("lame_lambda", "name of the lambda lame constant");
    params.addRequiredParam<MooseFunctorName>("lame_mu", "name of the mu lame constant");

    //also provide directly momentum components and density
    params.addRequiredParam<MooseFunctorName>("momentum_vector", "name of the momentum vector");
    params.addRequiredParam<MooseFunctorName>("density", "name of the density variable");

    ////
    params.addRequiredParam<MooseFunctorName>("stress_x_name", "name of the stress component in the x direction");
    params.addRequiredParam<MooseFunctorName>("stress_y_name", "name of the stress component in the y direction");
    params.addRequiredParam<MooseFunctorName>("stress_z_name", "name of the stress component in the z direction");

    params.addRequiredParam<MooseFunctorName>("T_x_name", "name of the mechanical flux tensor component in the x direction");
    params.addRequiredParam<MooseFunctorName>("T_y_name", "name of the mechanical flux tensor component in the y direction");
    params.addRequiredParam<MooseFunctorName>("T_z_name", "name of the mechanical flux tensor component in the z direction");
    return params;
    
}

ADComputeStressTensor::ADComputeStressTensor(const InputParameters &params)
    : FunctorMaterial(params),
      _grad_x(getFunctor<ADRealVectorValue>("map_x_gradient")),
      _grad_y(getFunctor<ADRealVectorValue>("map_y_gradient")),
      _grad_z(getFunctor<ADRealVectorValue>("map_z_gradient")),
      _lambda(getFunctor<ADReal>("lame_lambda")),
      _mu(getFunctor<ADReal>("lame_mu")),
      _m(getFunctor<ADRealVectorValue>("momentum_vector")),
      _rho(getFunctor<ADReal>("density"))

////////////////////////
{   
    //define a lambda that computes stress
    auto computeStress = [this](const auto & r, const auto & state) -> ADRankTwoTensor
    {
        //retrieve vector components
        const ADRealVectorValue grad_x = _grad_x(r, state);
        const ADRealVectorValue grad_y = _grad_y(r, state);
        const ADRealVectorValue grad_z = _grad_z(r, state);
        const ADReal lambda = _lambda(r, state);
        const ADReal mu = _mu(r, state);

        //initialize the tensor
        ADRankTwoTensor inv_F;
        inv_F.zero();

        //populate tensor
        for (unsigned int j = 0; j < 3; ++j){
            inv_F(0,j) = grad_x(j);
            inv_F(1,j) = grad_y(j);
            inv_F(2,j) = grad_z(j);
        }

        //this is the deformation gradient
        ADRankTwoTensor F = inv_F.inverse();
        ADReal J = F.det();

        //using this, we compute the required lagrangian deformation tensors
        ADRankTwoTensor C = F.transpose() * F;
        ADRankTwoTensor I;
        I.setToIdentity();

        //form stress
        ADRankTwoTensor PK2;
        PK2 = lambda * MetaPhysicL::log(J) * C.inverse() + mu * (I - C.inverse());

        //form sigma
        ADRankTwoTensor sigma = (1. / J) * (F * PK2 * F.transpose());
        return sigma;
    };
    
    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("stress_x_name"),
        [computeStress](const auto & r, const auto & state) -> ADRealVectorValue{
            //provide only component in x of stress tensor => vector
            const ADRankTwoTensor sigma_full = computeStress(r, state);

            //obtain only directional components
            ADRealVectorValue sigma_x;
            for (unsigned int j = 0; j < 3; ++j){
                sigma_x(j) = sigma_full(0,j);
            }
            return sigma_x;
            //
        }
    );

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("stress_y_name"),
        [computeStress](const auto & r, const auto & state) -> ADRealVectorValue{
            //provide only component in y of stress tensor => vector
            const ADRankTwoTensor sigma_full = computeStress(r, state);

            //obtain only directional components
            ADRealVectorValue sigma_y;
            for (unsigned int j = 0; j < 3; ++j){
                sigma_y(j) = sigma_full(1,j);
            }
            return sigma_y;
            //
        }
    );

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("stress_z_name"),
        [computeStress](const auto & r, const auto & state) -> ADRealVectorValue{
            //provide only component in z of stress tensor => vector
            const ADRankTwoTensor sigma_full = computeStress(r, state);

            //obtain only directional components
            ADRealVectorValue sigma_z;
            for (unsigned int j = 0; j < 3; ++j){
                sigma_z(j) = sigma_full(2,j);
            }
            return sigma_z;
            //
        }
    );

    //now we can compute the mechanical flux tensor
    auto computeReynoldsTensor = [this, computeStress](const auto & r, const auto & state) -> ADRankTwoTensor
    {
        //obtain velocity vector
        const ADRealVectorValue m = _m(r, state);
        const ADReal rho = _rho(r, state);

        //form initial term
        ADRankTwoTensor mcrossm;
        mcrossm.zero();

        for (unsigned int i = 0; i < 3; ++i){
            for (unsigned int j = 0; j < 3; ++j){
                mcrossm(i,j) = m(i) * m(j) / rho;
            }
        }

        //subtract stress
        const ADRankTwoTensor stress_full = computeStress(r, state);

        ADRankTwoTensor T = mcrossm - stress_full;
        return T;
    };

    //use the formed total tensor to expose the vector components for each direction

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("T_x_name"),
        [computeReynoldsTensor](const auto & r, const auto & state) -> ADRealVectorValue
        {
            //retrieve tensor
            const ADRankTwoTensor T_full = computeReynoldsTensor(r, state);

            //declare direction 
            ADRealVectorValue T_x;
            T_x.zero();

            //form x component
            for (unsigned int j = 0; j < 3; ++j){
                T_x(j) = T_full(0,j);
            }
            return T_x;
        }
    );

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("T_y_name"),
        [computeReynoldsTensor](const auto & r, const auto & state) -> ADRealVectorValue
        {
            //retrieve tensor
            const ADRankTwoTensor T_full = computeReynoldsTensor(r, state);

            //declare direction 
            ADRealVectorValue T_y;
            T_y.zero();

            //form x component
            for (unsigned int j = 0; j < 3; ++j){
                T_y(j) = T_full(1,j);
            }
            return T_y;
        }
    );

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("T_z_name"),
        [computeReynoldsTensor](const auto & r, const auto & state) -> ADRealVectorValue
        {
            //retrieve tensor
            const ADRankTwoTensor T_full = computeReynoldsTensor(r, state);

            //declare direction 
            ADRealVectorValue T_z;
            T_z.zero();

            //form x component
            for (unsigned int j = 0; j < 3; ++j){
                T_z(j) = T_full(2,j);
            }
            return T_z;
        }
    );
}