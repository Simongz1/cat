#pragma once

#include "ADTimeKernel.h"

// Forward Declarations
class ADNablaRhoVY;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class ADNablaRhoVY : public ADTimeKernel
{
public:
  ADNablaRhoVY(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

private:
  const ADMaterialProperty<Real> & _rho;
  const ADVariableValue & _vx;
  const ADVariableValue & _vy;
  //gradients
  const ADVariableGradient & _grad_vx;
  const ADVariableGradient & _grad_vy;
};