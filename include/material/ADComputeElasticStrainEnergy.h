#pragma once

#include "ADMaterial.h"
//#include "ElasticityTensorTools.h"
//#include "ADSingleVariableReturnMappingSolution.h"
//#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADComputeElasticStrainEnergy
  : public DerivativeMaterialInterface<ADMaterial>//,
    //public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADComputeElasticStrainEnergy(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const ADVariableValue &_temperature;
  const ADVariableValue &_c;
  const VariableName &_c_name;
  const ADMaterialProperty<RankTwoTensor> &_epsilon;
  const ADMaterialProperty<RankTwoTensor> &_epsilon_p;
  //const ADMaterialProperty<RankTwoTensor> &_epsilon_c;
  const ADMaterialProperty<RankTwoTensor> &_epsilon_T;
  const ADMaterialProperty<Real> &_alpha_thermal;
  const ADMaterialProperty<Real> &_density;
  const Real _lambda;
  const Real _mu;

  ADMaterialProperty<RankTwoTensor> &_epsilon_e;
  ADMaterialProperty<Real> &_We;
  ADMaterialProperty<RankTwoTensor> &_small_stress;
  const ADMaterialProperty<RankTwoTensor> &_F;
  const ADMaterialProperty<RankTwoTensor> &_C;
  ADMaterialProperty<RankTwoTensor> &_pk1_stress;
  ADMaterialProperty<RankTwoTensor> &_pk2_stress;
  ADMaterialProperty<RankTwoTensor> &_stress;

  ADMaterialProperty<Real> &_Hist;
  const MaterialProperty<Real> &_Hist_old;

  const ADMaterialProperty<Real> &_Wp;
  const ADMaterialProperty<Real> &_D;
  const std::string _D_name;
  const ADMaterialProperty<Real> &_dD;
  ADMaterialProperty<Real> &_dissipation;
  const ADMaterialProperty<RankTwoTensor> &_sp;

  const ADMaterialProperty<RankTwoTensor> &_epsilon_p_dot;
  //const ADMaterialProperty<RankTwoTensor> &_epsilon_c_dot;

  const Real _c1;
  const Real _c2;
  const Real _c3;

  ADMaterialProperty<Real> &_switch;

  ADMaterialProperty<Real> &_sM;
  ADMaterialProperty<RankTwoTensor> &_nMnM;
  ADMaterialProperty<Real> &_s_pressure;

  const Real _C0;
  const Real _C1;
  const ADMaterialProperty<Real> &_J_dot;
  const ADVariableValue &_h;
  const ADMaterialProperty<Real> &_sound_speed;
};