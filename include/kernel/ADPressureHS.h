#include "ADMatHeatSource.h"
#include <RankTwoTensor.h>
#include "RankFourTensor.h"
#include <vector>

//forward declarte the class object to be acted upon

class ADPressureHS : public ADMatHeatSource
{
public:
  ADPressureHS(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();


private:
  const VariableValue &_T;
  unsigned int _T_var_number;
  const MaterialProperty<RankFourTensor> &_elasticity_tensor;
  const VariableValue &_Yfinal;
  const ADMaterialProperty<Real> &_pressure_mg;
  const ADMaterialProperty<Real> &_pressure_JWL;
  const Real _beta_av;
  const ADMaterialProperty<Real> &_pressure_av;
  const ADMaterialProperty<Real> &_dP_dT;
  const MaterialProperty<RankTwoTensor> &_Ee_dot;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
};
