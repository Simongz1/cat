#include "ADComputeComponentVector.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeComponentVector);

InputParameters
ADComputeComponentVector::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("form vector from provided components");

    params.addRequiredParam<MooseFunctorName>("x_component", "x component");
    params.addRequiredParam<MooseFunctorName>("y_component", "y component");
    params.addRequiredParam<MooseFunctorName>("z_component", "z component");

    params.addRequiredParam<MooseFunctorName>("vector_name", "name of the vector field");
    return params;
}

ADComputeComponentVector::ADComputeComponentVector(const InputParameters &params)
    : FunctorMaterial(params),
      _x(getFunctor<ADReal>("x_component")),
      _y(getFunctor<ADReal>("y_component")),
      _z(getFunctor<ADReal>("z_component"))
{

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("vector_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{

            //forward declare the required quantities
            const ADReal x = _x(r, state);
            const ADReal y = _y(r, state);
            const ADReal z = _z(r, state);

            //form vector
            ADRealVectorValue out;
            out(0) = x;
            out(1) = y;
            out(2) = z;
            return out;
        }
    );
}