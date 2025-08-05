
#pragma once

#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class J2HeatSource : public HeatSource
{
public:
  
  static InputParameters validParams();
  J2HeatSource(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  // virtual Real computeQpJacobian();
  // virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  
  std::string _base_name;

  // reference temperature with zero thermal expansion
  const VariableValue & _temperature;
  const Real _reference_temperature;
  const MaterialProperty<Real> & _specific_heat;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor; //elasticity tensor
  const Real _alpha;
  const Real _eta;
  const MaterialProperty<Real> & _effective_plastic_strain_old;
  const MaterialProperty<Real> & _effective_plastic_strain;
  const VariableValue & _effective_ps;
  const VariableValue & _vol_strain;
  const MaterialProperty<RankTwoTensor> & _pk1_stress;
  // const MaterialProperty<RankTwoTensor> & _pk2_stress;
  const Real _beta_p;
  const MaterialProperty<RankTwoTensor> & _total_strain_old;
  const MaterialProperty<RankTwoTensor> & _total_strain;
  const Real _gamma; // gruneisen parameter
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  const Real _C0;
  const Real _C1;
  const Real _h;
  const Real _ss;
  const Real _beta_vis;
  
  // MaterialProperty<Real> & _effective_ps,

};

