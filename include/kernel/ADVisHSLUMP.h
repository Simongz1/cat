#include "ADMatHeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "ElasticityTensorTools.h"
#include <vector>

//forward declarte the class object to be acted upon

class ADVisHSLUMP : public ADMatHeatSource
{
public:
  ADVisHSLUMP(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();
  

private:

  const MaterialProperty<RankTwoTensor> &_cauchy_stress;
  const Real _beta_p;
  const MaterialProperty<RankTwoTensor> &_Ep_dot;
  const ADMaterialProperty<Real> &_rho;
  const ADMaterialProperty<Real> &_cv;
  const MaterialProperty<RankTwoTensor> &_F;
  const Real _use_PK2;
  const VariableValue _dirac_switch_react;
};
