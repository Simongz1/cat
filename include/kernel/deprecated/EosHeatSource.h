#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class EosHeatSource;

//template <>
//InputParameters validParams<ThermalExpansionHeatSourceFiniteStrainMieGruneisenNew>();

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)
 */
class EosHeatSource : public HeatSource
{
public:
  EosHeatSource(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  std::string _base_name;

  // Gruneisen G (or Gamma) parameter (Menon, 2014)
  const Real _G_Gruneisen;

  // Reference bulk modulus in the equation of state
  const Real _Bulk_Modulus_Ref;

  // reference temperature with zero thermal expansion
  const Real _reference_temperature;

  const MaterialProperty<Real> & _specific_heat;
  const MaterialProperty<Real> & _density;
  

  const MaterialProperty<RankTwoTensor> & _deformation_gradient; // deformation gradient
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old; // deformation gradient, previous timestep

  const MaterialProperty<RankFourTensor> & _elasticity_tensor; //elasticity tensor
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  const MaterialProperty<RankTwoTensor> & _total_strain;
  const MaterialProperty<RankTwoTensor> & _total_strain_old;
  const MaterialProperty<RankTwoTensor> & _elastic_strain;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;

  // MaterialProperty<RankTwoTensor> & _PK1;
  // MaterialProperty<RankTwoTensor> & _PK2;

  const MaterialProperty<RankTwoTensor> & _stress;
  const Real _beta_av;
  const Real _beta_p;

  const Real _C0;
  const Real _C1;
  const Real _h_max;

  const Real _ss;

  //MaterialProperty<RankTwoTensor> & _PK1;
  //MaterialProperty<RankTwoTensor> & _PK2;

  //MaterialProperty<Real> & _q_eos;
  //MaterialProperty<Real> & _q_cpl;
  //MaterialProperty<Real> & _q_av;
  //MaterialProperty<Real> & _q_p;
  //MaterialProperty<Real> & _q_tot;
};

