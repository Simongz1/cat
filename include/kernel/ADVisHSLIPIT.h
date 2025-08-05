#include "ADMatHeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "ElasticityTensorTools.h"
#include <vector>

//forward declarte the class object to be acted upon

class ADVisHSLIPIT : public ADMatHeatSource
{
public:
  ADVisHSLIPIT(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();
  

private:

  const MaterialProperty<RankTwoTensor> &_cauchy_stress;
  const Real _beta_p;
  const Real _beta_comp;
  const MaterialProperty<RankTwoTensor> &_Ep_dot;
  const MaterialProperty<RankTwoTensor> &_Ee_dot;
  const MaterialProperty<RankTwoTensor> &_sigma;
  const MaterialProperty<RankTwoTensor> &_F;
  const MaterialProperty<Real> &_ep_rate;
  const ADMaterialProperty<Real> &_alpha;
  const MaterialProperty<RankTwoTensor> &_Fe;
  const MaterialProperty<RankTwoTensor> &_Fe_old;
  const Real _Tref;
  const MaterialProperty<RankFourTensor> &_Cijkl;
};
