#include "FVInterpolatedDivergence.h"
#include "FVFluxKernel.h"

registerADMooseObject("mlApp", FVInterpolatedDivergence);

InputParameters
FVInterpolatedDivergence::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params.addClassDescription("enables tunable selection of interpolation method for divergence enforcement.");
  params.addParam<std::string>("interpolation_method", "Average", "interpolation method to use for face values");
  params.addRequiredParam<MooseFunctorName>("vector_field", "name of the vector field to apply divergence to");
  return params;
}

FVInterpolatedDivergence::FVInterpolatedDivergence(const InputParameters & parameters)
  : FVFluxKernel(parameters),
    _vector_field(getFunctor<ADRealVectorValue>("vector_field")),
    _interpolation_method(getParam<std::string>("interpolation_method"))
{}

ADReal
FVInterpolatedDivergence::computeQpResidual()
{ 
  //make face with the custom interpolation method
  Moose::FV::InterpMethod method;
  //for now only support average or upwind
  if (_interpolation_method != "Average"){
     method = Moose::FV::InterpMethod::Upwind;
  }
  else{
    method = Moose::FV::InterpMethod::Average;
  }

  const auto face = makeFace(*_face_info, Moose::FV::limiterType(method), true);

  //get vector
  const ADRealVectorValue vector = _vector_field(face, determineState());

  //return statement
  return 1.0 * (vector * _normal);
}