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

class ComputeMieGruneisenPlasticPFFractureStress : public ComputeStressBase
{
public:
  static InputParameters validParams();
  ComputeMieGruneisenPlasticPFFractureStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress() override;
  virtual void initQpStatefulProperties() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  const VariableValue & _c;
  const MaterialProperty<Real> & _l;
  const MaterialProperty<Real> & _pressure;
  const MaterialProperty<Real> & _gc;

  bool _use_current_hist;
  bool _use_snes_vi_solver;

  MaterialProperty<Real> & _H;
  const MaterialProperty<Real> & _H_old;
  const MaterialProperty<Real> & _barrier;
  MaterialProperty<Real> & _E;
  MaterialProperty<Real> & _dEdc;
  MaterialProperty<Real> & _d2Ed2c;
  MaterialProperty<RankTwoTensor> & _dstress_dc;
  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;
  const MaterialProperty<Real> & _D;
  const MaterialProperty<Real> & _dDdc;
  const MaterialProperty<Real> & _d2Dd2c;
  const MaterialProperty<Real> & _I;
  const MaterialProperty<Real> & _dIdc;
  const MaterialProperty<Real> & _d2Id2c;
  const Real _Gamma;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<Real> & _specific_heat;
  const VariableValue & _temperature;
  const Real _ref_temperature;
  const Real _s;

  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  const Real _C0;
  const Real _C1;
  const Real _Le;
  const Real _sound_speed;

  MaterialProperty<Real> & _bulk_modulus;
  MaterialProperty<Real> & _pressure_eos;

  std::vector<Real> _yield_stress_vector;
  MaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
  MaterialProperty<Real> & _eqv_plastic_strain;
  const MaterialProperty<Real> & _eqv_plastic_strain_old;
  const MaterialProperty<RankTwoTensor> & _rotation_increment;
  Real _rtol;
  Real _ftol;
  Real _eptol;

  // outer and mixed product of the delta function tensors
  RankFourTensor _deltaOuter, _deltaMixed;
  MaterialProperty<Real> & _W0p;

  MaterialProperty<RankTwoTensor> & _stress_eos_elastic;
  MaterialProperty<RankTwoTensor> & _stress_cpl_elastic;
  MaterialProperty<RankTwoTensor> & _stress_eos_plastic;
  MaterialProperty<RankTwoTensor> & _stress_cpl_plastic;

  virtual Real yieldFunction(const RankTwoTensor & stress, const Real yield_stress);
  virtual RankTwoTensor dyieldFunction_dstress(const RankTwoTensor & stress);
  virtual Real dyieldFunction_dinternal(const Real equivalent_plastic_strain);
  virtual RankTwoTensor flowPotential(const RankTwoTensor & stress);
  virtual Real internalPotential();
  Real getSigEqv(const RankTwoTensor & stress);
  virtual void getJac(const RankTwoTensor & sig,
                      const RankFourTensor & E_ijkl,
                      Real flow_incr,
                      RankFourTensor & dresid_dsig);
  Real getYieldStress(const Real equivalent_plastic_strain);
  Real getdYieldStressdPlasticStrain(const Real equivalent_plastic_strain);
  usingTensorIndices(i_, j_, k_, l_);
};
