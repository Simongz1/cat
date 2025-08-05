/// Calculates heat generated due to thermal expansion

#pragma once

#include "HeatSource.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Luscher2017 model and Amor damage model
 * Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)
 */
class ThermalExpansionHeatSourceSmallStrainMieGruneisen : public HeatSource
{
public:
  static InputParameters validParams();

  ThermalExpansionHeatSourceSmallStrainMieGruneisen(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  std::string _base_name;

  // Gruneisen G (or Gamma) parameter (Menon, 2014)
  const Real _Gamma;

  const MaterialProperty<Real> & _specific_heat;
  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  // damage phase field
  const VariableValue & _c;

  const bool _c_coupled;
  const unsigned int _c_var;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

  const MaterialProperty<Real> & _shock_heat;
  // const MaterialProperty<Real> & _expansion_heat;
  // const MaterialProperty<Real> & _vis_heat;

};
