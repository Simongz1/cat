/// Calculates heat generated due to thermal expansion

#pragma once

#include "Material.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class ADComputeMISTERnetHeat : public Material
{
public:
  static InputParameters validParams();

  ADComputeMISTERnetHeat(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  std::string _base_name;
  const Real _T_ref;
  const VariableValue & _dirac_switch_shock;
  const VariableValue & _dirac_switch_react;
  const VariableValue & _T;

  const ADMaterialProperty<Real> & _density;
  const ADMaterialProperty<Real> & _specific_heat;

  const Real _factor;
  const Real _heat_time_shock;
  const Real _heat_time_react;
  const MaterialProperty<Real> & _temperature_mister_shock;
  const MaterialProperty<Real> & _temperature_mister_react;
  ADMaterialProperty<Real> & _heatrate_mister_shock;
  ADMaterialProperty<Real> & _heatrate_mister_react;

  const MaterialProperty<Real> &_v_flag;
  const VariableValue &_vx;
  const VariableValue &_ax;
  const Real _thr_v;
  const Real _thr_a;
};
