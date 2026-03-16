#pragma once

#include "ADMaterial.h"
//#include "ElasticityTensorTools.h"
//#include "ADSingleVariableReturnMappingSolution.h"
//#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADComputePlasticFlowRate
  : public DerivativeMaterialInterface<ADMaterial>//,
    //public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADComputePlasticFlowRate(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const ADVariableValue &_temperature;
  const ADVariableValue &_shear_strength;
  const Real _gamma_p_0;
  const Real _alpha_shear;
  const Real _h;
  const Real _ss;

  const Real _A;
  const ADMaterialProperty<RankTwoTensor> &_sp;
  const ADMaterialProperty<Real> &_s_pressure;

  ADMaterialProperty<Real> &_gamma_p_dot;
  ADMaterialProperty<Real> &_shear_strength_dot;
};