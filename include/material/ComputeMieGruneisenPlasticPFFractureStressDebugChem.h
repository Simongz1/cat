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

class ComputeMieGruneisenPlasticPFFractureStressDebugChem : public ComputeStressBase
{
public:
  static InputParameters validParams();
  ComputeMieGruneisenPlasticPFFractureStressDebugChem(const InputParameters & parameters);

protected:
  virtual void computeQpStress() override;
  virtual void initQpStatefulProperties() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

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
  // const Real _sound_speed;

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
   
  MaterialProperty<Real> & _eqv_trial;
  const Real _yield3_thr;
  MaterialProperty<Real> & _dyield1;
  MaterialProperty<Real> & _ss_prop;
  const Real _gamma_pd;

  //chemical reaction parameters

  const VariableValue & _chemY1;
  const VariableValue & _chemY2;
  const VariableValue & _chemY3;
  const VariableValue & _chemY4;
  MaterialProperty <Real> & _chemY1_dot;
  MaterialProperty <Real> & _chemY2_dot;
  MaterialProperty <Real> & _chemY3_dot;
  MaterialProperty <Real> & _chemY4_dot;

  const Real _Z1;
  const Real _Z2;
  const Real _Z3;
  const Real _Z4;

  const Real _E1;
  const Real _E2;
  const Real _E3;

  const Real _R_const;

  const Real _A_g;
  const Real _B_g;
  const Real _C_g;
  const Real _R_1g;
  const Real _R_2g;
  const Real _omega;

  const Real _A_s;
  const Real _B_s;
  const Real _C_s;
  const Real _R_1s;
  const Real _R_2s;
  const Real _omega_s;

  MaterialProperty<Real> & _p_g;
  MaterialProperty<Real> & _p_s;
  
  //MaterialProperty<Real> & _lambda_prop;
  MaterialProperty<Real> & _lambda_dot;
  MaterialProperty<Real> & _lambda;
  const MaterialProperty<Real> & _lambda_old;
  //const MaterialProperty<Real> & _lambda_prop_old;
  //const VariableValue & _lambda; //variable definition of lambda

  const Real _I_ev;
  const Real _G1;
  const Real _G2;

  const Real _x0;
  const Real _y0;
  const Real _x1;
  const Real _y1;
  const Real _z1;
  const Real _x2;
  const Real _y2;
  const Real _z2;
  const Real _a0;

  const Real _useArrhenius;

  MaterialProperty<Real> & _density_corr;
  MaterialProperty<Real> & _mu;

  const Real _isTempDep;
  const Real _isMG;

  const MaterialProperty<Real> & _Cv_gas;
  const MaterialProperty<Real> & _Cv_solid;

  MaterialProperty<Real> & _lambda1;
  MaterialProperty<Real> & _lambda2;
  MaterialProperty<Real> & _lambda3;

  const Real _igmax;
  const Real _gmax1;

  MaterialProperty<Real> & _pressure_eos_plastic;
  MaterialProperty<Real> & _jac_lambda;

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
  // virtual Real chem_react_1(const Real Y1, const Real T);
  // virtual Real chem_react_2(const Real Y1, const Real Y1_rate, const Real Y2, const Real T);
  // virtual Real chem_react_3(const Real Y2, const Real Y2_rate, const Real Y3, const Real T);
  // virtual Real chem_react_4(const Real Y3, const Real Y3_rate, const Real T);
  usingTensorIndices(i_, j_, k_, l_);
};
