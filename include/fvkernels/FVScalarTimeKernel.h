#pragma once

#include "FVTimeKernel.h"

class FVScalarTimeKernel : public FVTimeKernel
{
public:
  static InputParameters validParams();
  FVScalarTimeKernel(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;
  const Moose::Functor<ADReal> & _scalar;
};