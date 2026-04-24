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
    return params;
}

ADComputeUVector::ADComputeUVector(const InputParameters &params)
    : FunctorMaterial(params),
      _density(getFunctor<ADReal>("density")),
      _m_vector(getFunctor<ADRealVectorValue>("m_vector"))
{

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("u_vector_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{

            //forward declare the required quantities
            const ADReal density = _density(r, state);
            const ADRealVectorValue m_vector = _m_vector(r, state);

            //define the output
            return m_vector / density; //velocity vector
        }
    );
}