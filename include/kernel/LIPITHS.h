#include "ADKernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "ElasticityTensorTools.h"
#include <vector>

//forward declarte the class object to be acted upon

class LIPITHS : public ADKernel
{
public:
  LIPITHS(const InputParameters & parameters);
  static InputParameters validParams();
protected:
  virtual ADReal computeQpResidual() override;
private:
  const Real _beta_p;
  const Real _beta_av;
  const ADMaterialProperty<RankTwoTensor> &_Ep_dot;
  const ADMaterialProperty<RankTwoTensor> &_Ee_dot;
  const ADMaterialProperty<Real> &_alpha;
  const ADMaterialProperty<RankTwoTensor> &_F;

  const ADMaterialProperty<RankFourTensor> &_Cijkl;
  const ADMaterialProperty<RankTwoTensor> &_S;
  const ADMaterialProperty<Real> &_HS_plastic;
  const ADMaterialProperty<Real> &_HS_elastic;
};

