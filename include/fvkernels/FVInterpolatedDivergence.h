#pragma once

#include "FVFluxKernel.h"

class FVInterpolatedDivergence : public FVFluxKernel
{
public:
  static InputParameters validParams();
  FVInterpolatedDivergence(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;
  const Moose::Functor<ADRealVectorValue> &_vector_field;
  const std::string _interpolation_method;
};