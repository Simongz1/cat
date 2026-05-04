#include "ADComputeVectorScalar.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeVectorScalar);

InputParameters
ADComputeVectorScalar::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes a vector field as the result of an input vector field times a scalar functor");

    params.addRequiredParam<MooseFunctorName>("scalar", "name of the scalar quantity");
    params.addRequiredParam<MooseFunctorName>("vector", "name of the vector quantity");
    params.addRequiredParam<MooseFunctorName>("u_vector_name", "name of the U vector field");
    return params;
}

ADComputeVectorScalar::ADComputeVectorScalar(const InputParameters &params)
    : FunctorMaterial(params),
      _scalar(getFunctor<ADReal>("scalar")),
      _vector(getFunctor<ADRealVectorValue>("vector"))
{

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("u_vector_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{

            //forward declare the required quantities
            const ADReal scalar = _scalar(r, state);
            const ADRealVectorValue vector = _vector(r, state);

            //define the output
            return scalar * vector; //velocity vector
        }
    );
}