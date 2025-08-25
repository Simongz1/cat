/// Calculates heat generated due to thermal expansion

#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class ComputeMolecularJettingHeatSource : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ComputeMolecularJettingHeatSource(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  std::string _base_name;

  // damage phase field
  const VariableValue & _c;
  const VariableValue & _dirac_switch;

  const MaterialProperty<Real> & _density;
  const MaterialProperty<Real> & _specific_heat;

  const Real _alpha;
  const Real _up;
  const Real _s;
  const Real _sound_speed;
  const Real _g;
  const Real _zeta;
  const Real _tau_induction;
  const Real _factor;
  MaterialProperty<Real> & _molecular_jetting_heat;
  MaterialProperty<Real> & _d_molecular_jetting_heat_dT;
  MaterialProperty<Real> & _d_molecular_jetting_heat_dc;
};
