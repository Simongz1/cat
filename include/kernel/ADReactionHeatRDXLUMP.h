#include "ADMatHeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <vector>

//forward declarte the class object to be acted upon

class ADReactionHeatRDXLUMP : public ADMatHeatSource
{
public:
  ADReactionHeatRDXLUMP(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();

private:
  const ADMaterialProperty<Real> &_Q1;
  const ADMaterialProperty<Real> &_Q2;
  const ADMaterialProperty<Real> &_r1;
  const ADMaterialProperty<Real> &_r2;

  const VariableValue & _Y1;
  const VariableValue & _Y2;
  const VariableValue & _Y3;
  const ADMaterialProperty<Real> &_rho;
  const ADMaterialProperty<Real> &_cv;
  const VariableValue &_dirac_switch_react;
};
