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

class ComputeMieGruneisenPlasticPFFractureStressNew : public ComputeStressBase
{
public:
  static InputParameters validParams();
  ComputeMieGruneisenPlasticPFFractureStressNew(const InputParameters & parameters);

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

  // std::vector<Real> _yield_stress_vector;
  MaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
  MaterialProperty<RankTwoTensor> & _plastic_strain_rate;
  MaterialProperty<Real> & _eqv_plastic_strain;
  MaterialProperty<Real> & _eqv_plastic_strain_rate;
  const MaterialProperty<Real> & _eqv_plastic_strain_rate_old;
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

  const Real _A_yield;
  const Real _B_yield;
  const Real _C_yield;
  const Real _a_corr_yield;
  const Real _m_yield;
  const Real _n_yield;
  const Real _strain_ref;
  const Real _T_melt_ref;
  const Real _maxiter;
  const Real _toly;
  const Real _k_alpha;
  const Real _thr_gas;
  
  MaterialProperty<Real> & _yield_1;
  MaterialProperty<Real> & _yield_2;
  MaterialProperty<Real> & _yield_3;
  MaterialProperty<Real> & _yield_total;
  MaterialProperty<Real> & _dyield_dplastic;
  MaterialProperty<Real> & _theta_norm;
  MaterialProperty<Real> & _T_melt;
  MaterialProperty<Real> & _V;
  
  const Real _thr_linear;

  virtual Real yieldFunction(const RankTwoTensor & stress, const Real yield_stress);
  virtual RankTwoTensor dyieldFunction_dstress(const RankTwoTensor & stress);
  virtual Real dyieldFunction_dinternal(const Real eqv_plastic_strain, const Real eqv_plastic_strain_rate, const Real Temperature);
  virtual RankTwoTensor flowPotential(const RankTwoTensor & stress);
  virtual Real internalPotential();
  Real getSigEqv(const RankTwoTensor & stress);
  virtual void getJac(const RankTwoTensor & sig,
                      const RankFourTensor & E_ijkl,
                      Real flow_incr,
                      RankFourTensor & dresid_dsig);
  Real getYieldStress(const Real eqv_plastic_strain, const Real eqv_plastic_strain_rate, const Real Temperature);
  Real getdYieldStressdPlasticStrain(const Real eqv_plastic_strain, const Real eqv_plastic_strain_rate, const Real Temperature);
  usingTensorIndices(i_, j_, k_, l_);
};
