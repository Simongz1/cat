#pragma once

#include "ADMaterial.h"
//#include "ElasticityTensorTools.h"
//#include "ADSingleVariableReturnMappingSolution.h"
//#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADComputeCrazingFlowRate
  : public DerivativeMaterialInterface<ADMaterial>//,
    //public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADComputeCrazingFlowRate(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const Real _gamma_c_0;
  const Real _m_craze;
  const Real _craze_flow_resistance;
  const ADMaterialProperty<Real> &_sM;

  ADMaterialProperty<Real> &_gamma_c_dot;
};