#pragma once

#include "ADKernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <vector>

//forward declarte the class object to be acted upon

class ADQYdot : public ADKernel
{
public:
  ADQYdot(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();

private:
  const Real &_Q1;
  const Real &_Q2;
  const Real &_Q3;

  const ADMaterialProperty<Real> &_Y1dot;
  const ADMaterialProperty<Real> &_Y2dot;
  const ADMaterialProperty<Real> &_Y3dot;
  const ADMaterialProperty<Real> &_Y4dot;
};