#pragma once

#include "ADMaterial.h"
//#include "ElasticityTensorTools.h"
//#include "ADSingleVariableReturnMappingSolution.h"
//#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADIntegrateStrains
  : public DerivativeMaterialInterface<ADMaterial>//,
    //public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADIntegrateStrains(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  ADMaterialProperty<RankTwoTensor> &_epsilon_p;
  //ADMaterialProperty<RankTwoTensor> &_epsilon_c;
  ADMaterialProperty<RankTwoTensor> &_epsilon_T;

  const MaterialProperty<RankTwoTensor> &_epsilon_p_old;
  //const MaterialProperty<RankTwoTensor> &_epsilon_c_old;

  const MaterialProperty<RankTwoTensor> &_epsilon_p_dot_old;
  //const ADMaterialProperty<RankTwoTensor> &_epsilon_c_dot;

  const ADVariableValue &_temperature;
  const ADMaterialProperty<Real> &_alpha_thermal;
};