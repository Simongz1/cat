/// Calculates heat generated due to thermal expansion

#pragma once

#include "Material.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <cmath>

class ADComputeMISTERnetHeat : public Material
{
public:
  static InputParameters validParams();

  ADComputeMISTERnetHeat(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual ADReal CosGate(const ADReal tau, const ADReal t);

  std::string _base_name;
  const Real _T_ref;
  const ADVariableValue & _dirac_switch_shock;
  const ADVariableValue & _dirac_switch_react;
  const ADVariableValue & _T;

  const ADMaterialProperty<Real> & _density;
  const ADMaterialProperty<Real> & _specific_heat;

  const MaterialProperty<Real> & _temperature_mister_shock;
  const MaterialProperty<Real> & _temperature_mister_react;
  ADMaterialProperty<Real> & _heatrate_mister_shock;
  ADMaterialProperty<Real> & _heatrate_mister_react;

  const ADMaterialProperty<Real> &_v_flag;
  const std::vector<VariableName> _v_components;
  const std::vector<VariableName> _a_components;

  //surrogate chemistry source

  ADMaterialProperty<Real> &_Y1_dot_surrogate;
  ADMaterialProperty<Real> &_Y2_dot_surrogate;
  ADMaterialProperty<Real> &_Y3_dot_surrogate;
  MaterialProperty<Real> &_indicator_surrogate;
  const ADMaterialProperty<Real> &_time_react;
  ADMaterialProperty<Real> &_time_shock;
  const Real _h;
  const VariableValue & _fraction_csv;

  const bool _correction_heat;
  const ADMaterialProperty<Real> &_us;
  const bool _use_gating;
  const bool _use_complete_burn;

  std::vector<const ADVariableValue *> _v;
  std::vector<const ADVariableValue *> _a;
};
