#pragma once

#include "ADMaterial.h"
//#include "ElasticityTensorTools.h"
//#include "ADSingleVariableReturnMappingSolution.h"
//#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADComputeCrazingStrainRate
  : public DerivativeMaterialInterface<ADMaterial>//,
    //public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADComputeCrazingStrainRate(const InputParameters & parameters);

  virtual void initialSetup() override;
  

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  const ADMaterialProperty<Real> &_gamma_c_dot;
  const ADMaterialProperty<Real> &_switch;
  const ADMaterialProperty<RankTwoTensor> &_nMnM;
  ADMaterialProperty<RankTwoTensor> &_epsilon_c_dot;
};