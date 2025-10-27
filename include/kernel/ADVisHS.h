#include "ADMatHeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "ElasticityTensorTools.h"
#include <vector>

//forward declarte the class object to be acted upon

class ADVisHS : public ADMatHeatSource
{
public:
  ADVisHS(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();
  

private:

  const MaterialProperty<RankTwoTensor> &_cauchy_stress;
  const Real _beta_p;
  const MaterialProperty<RankTwoTensor> &_Ep_dot;
};
