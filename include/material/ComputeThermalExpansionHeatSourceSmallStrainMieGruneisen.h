/// Calculates heat generated due to thermal expansion

#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Luscher2017 model and Amor damage model
 * Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)
 */
class ComputeThermalExpansionHeatSourceSmallStrainMieGruneisen : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ComputeThermalExpansionHeatSourceSmallStrainMieGruneisen(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  std::string _base_name;

  /// Name of the elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  // Gruneisen G (or Gamma) parameter (Menon, 2014)
  const Real _Gamma;

  // reference temperature with zero thermal expansion
  const Real _ref_temperature;

  const MaterialProperty<Real> & _specific_heat;
  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  // damage phase field
  const VariableValue & _c;
  const VariableValue & _temp;

  /// Maximum element size
  const Real _Le;

  const bool _c_coupled;
  const unsigned int _c_var;

  // Von Neumann coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

  /// Material property defining speed of sound in material
  const Real _sound_speed;

  const Real _beta_v;

  MaterialProperty<Real> & _shock_heat;
  MaterialProperty<Real> & _expansion_heat;
  MaterialProperty<Real> & _expansion_heat_large;
  MaterialProperty<Real> & _vis_heat;

  MaterialProperty<RankTwoTensor> & _invce;
  MaterialProperty<RankTwoTensor> & _ce;
  MaterialProperty<Real> & _jacob;
  MaterialProperty<RankTwoTensor> & _ee_rate;
  MaterialProperty<RankTwoTensor> & _invce_ee_rate;

 

};
