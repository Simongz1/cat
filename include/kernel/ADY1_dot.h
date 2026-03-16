#pragma once

#include "ADKernel.h"

// Forward Declarations
class ADY1_dot;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class ADY1_dot : public ADKernel
{
public:
  ADY1_dot(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual() override;

private:
  std::string _base_name;
  const ADMaterialProperty<Real> &_Y1dot;
};

