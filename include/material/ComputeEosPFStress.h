/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State

#pragma once

#include "ComputeStressBase.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "PiecewiseBilinear.h"
#include "RankFourTensor.h"
#include "ExternalPetscSolverApp.h"
#include "petscblaslapack.h"

class ComputeEosPFStress : public ComputeStressBase
{
public:
  static InputParameters validParams();
  ComputeEosPFStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress() override;
  virtual void initQpStatefulProperties() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  // const MaterialProperty<Real> & _pressure;
 
  const Real _Gamma;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<Real> & _specific_heat;

// temperature
  const VariableValue & _temperature;
  const Real _ref_temperature;


  const Real _s;

  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
  const MaterialProperty<RankTwoTensor> & _total_strain;
  const MaterialProperty<RankTwoTensor> & _total_strain_old;

  const Real _C0;
  const Real _C1;
  const Real _Le;
  const Real _ss;

  MaterialProperty<Real> & _bulk_modulus;
  MaterialProperty<Real> & _pressure_eos;
  MaterialProperty<Real> & _peos1;
  MaterialProperty<Real> & _peos2;

  const VariableValue & _c;
  const Real _k_dmg;
  MaterialProperty<Real> & _E;
  MaterialProperty<Real> & _dEdc;
  MaterialProperty<Real> & _d2Ed2c;
  const MaterialProperty<Real> & _I;
  const MaterialProperty<Real> & _dIdc;
  const MaterialProperty<Real> & _d2Id2c;
  const MaterialProperty<Real> & _D;
  const MaterialProperty<Real> & _dDdc;
  const MaterialProperty<Real> & _d2Dd2c;

  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;

  const MaterialProperty<Real> & _f_pressure;
  MaterialProperty<Real> & _H;
  const MaterialProperty<Real> & _H_old;
  bool _use_current_hist;
  MaterialProperty<Real> & _elastic_energy;
  MaterialProperty<RankTwoTensor> & _deviatoric_strain;
  // MaterialProperty<Real> & _elastic_energy;





  usingTensorIndices(i_, j_, k_, l_);
};