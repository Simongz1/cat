#pragma once

#include "FVDiffusion.h"

class FVCoupledDiffusion : public FVDiffusion
{
public:
  static InputParameters validParams();
  FVCoupledDiffusion(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;
  const MooseVariableFV<Real> &_coupled_variable;
  const Real _factor;
};