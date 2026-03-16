#pragma once

#include "ADMaterial.h"
//#include "ElasticityTensorTools.h"
//#include "ADSingleVariableReturnMappingSolution.h"
//#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADComputePlasticStrainRate
  : public DerivativeMaterialInterface<ADMaterial>//,
    //public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADComputePlasticStrainRate(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const ADMaterialProperty<Real> &_gamma_p_dot;
  const ADMaterialProperty<RankTwoTensor> &_sp;
  const ADMaterialProperty<Real> &_switch;

  ADMaterialProperty<RankTwoTensor> &_epsilon_p_dot;
};