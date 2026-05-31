#include "ADComputeMaterialPointLocation.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>
#include "MathFVUtils.h"
#include <type_traits>
#include "MooseFunctorArguments.h"
#include "FaceInfo.h"

registerMooseObject("mlApp", ADComputeMaterialPointLocation);

InputParameters
ADComputeMaterialPointLocation::validParams(){

    InputParameters params = FunctorMaterial::validParams();
    params.addClassDescription("reconstructs the Lagrangian deformation gradient from the partial gradients of each component of the inverse map field");

    //request the full material point position vector
    params.addRequiredParam<MooseFunctorName>("rhoX_name", "name of the material point position vector multiplied by the density");
    params.addRequiredParam<MooseFunctorName>("X_name", "name of the material point position vector without density");
    params.addRequiredParam<MooseFunctorName>("density_name","name of the density variable");
    return params;
    
}

ADComputeMaterialPointLocation::ADComputeMaterialPointLocation(const InputParameters &params)
    : FunctorMaterial(params),
      _rhoX(getFunctor<ADRealVectorValue>("rhoX_name")),
      _rho(getFunctor<ADReal>("density_name"))

////////////////////////
{   
    //declare the material property needed
    addFunctorProperty<ADRealVectorValue>(getParam<MooseFunctorName>("X_name"),
        [this](const auto r, const auto state) -> ADRealVectorValue{
            return _rhoX(r, state) / _rho(r, state);
        });

    //use the same object to obtain the the gradient of the material point
    
}