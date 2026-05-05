#include "ADComputeUVector.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeUVector);

InputParameters
ADComputeUVector::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes the vector u * J for artificial viscosity advection");

    params.addRequiredParam<MooseFunctorName>("density", "name of the density variable");
    params.addRequiredParam<MooseFunctorName>("m_vector", "momentum vector");
    params.addRequiredParam<MooseFunctorName>("u_vector_name", "name of the U vector field");

    //add a sign parameter to switch from positive to negative
    params.addParam<Real>("scalar_sign", 1., "sign scalar for switching from positive to negative vector definition");
    return params;
}

ADComputeUVector::ADComputeUVector(const InputParameters &params)
    : FunctorMaterial(params),
      _density(getFunctor<ADReal>("density")),
      _m_vector(getFunctor<ADRealVectorValue>("m_vector")),
      _sign(getParam<Real>("scalar_sign"))
{

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("u_vector_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{

            //forward declare the required quantities
            const ADReal density = _density(r, state);
            const ADRealVectorValue m_vector = _m_vector(r, state);
            const Real sign = _sign;

            //define the output
            return sign * m_vector / density; //velocity vector
        }
    );
}