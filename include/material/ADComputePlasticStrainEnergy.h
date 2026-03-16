#pragma once

#include "ADMaterial.h"
//#include "ElasticityTensorTools.h"
//#include "ADSingleVariableReturnMappingSolution.h"
//#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADComputePlasticStrainEnergy
  : public DerivativeMaterialInterface<ADMaterial>//,
    //public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADComputePlasticStrainEnergy(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  virtual ADReal computeInverseLang(const ADReal &x);
  virtual ADReal computeInverseLangDerivative(const ADReal &x);
  virtual ADReal computedWdInvLang(const ADReal &x);
  virtual ADRankTwoTensor computedLambdadep();


  const ADMaterialProperty<RankTwoTensor> &_epsilon_p;
  const Real _N_p;
  const Real _mu_p;

  ADMaterialProperty<RankTwoTensor> &_exp2ep;
  ADMaterialProperty<Real> &_Wp;
  ADMaterialProperty<RankTwoTensor> &_sp;
};