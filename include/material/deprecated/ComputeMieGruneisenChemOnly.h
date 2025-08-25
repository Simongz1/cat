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

class ComputeMieGruneisenChemOnly : public ComputeStressBase
{
public:
  static InputParameters validParams();
  ComputeMieGruneisenChemOnly(const InputParameters & parameters);

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

  MaterialProperty<Real> & _bulk_modulus;
  MaterialProperty<Real> & _ss_prop;

  //chemical reaction parameters

  MaterialProperty<Real> & _chemY1;
  MaterialProperty<Real> & _chemY2;
  MaterialProperty<Real> & _chemY3;
  MaterialProperty<Real> & _chemY4;
  
  MaterialProperty<Real> & _chemY1_dot;
  MaterialProperty<Real> & _chemY2_dot;
  MaterialProperty<Real> & _chemY3_dot;
  MaterialProperty<Real> & _chemY4_dot;

  const Real _Z1;
  const Real _Z2;
  const Real _Z3;
  const Real _Z4;

  const Real _E1;
  const Real _E2;
  const Real _E3;

  const Real _R_const;

  const MaterialProperty<Real> & _chemY1_old;
  const MaterialProperty<Real> & _chemY2_old;
  const MaterialProperty<Real> & _chemY3_old;
  const MaterialProperty<Real> & _chemY4_old;

  const Real _A_g;
  const Real _B_g;
  const Real _C_g;
  const Real _R_1g;
  const Real _R_2g;
  const Real _omega;

  MaterialProperty<Real> & _p_g;
  MaterialProperty<Real> & _p_s;
  MaterialProperty<Real> & _lambda_dot;
  MaterialProperty<Real> & _lambda;

  const MaterialProperty<Real> & _lambda_old;

  const Real _I_ev;
  const Real _G_ev;
  const Real _a_ev;
  const Real _b_ev;

  const Real _useArrhenius;

  MaterialProperty<Real> & _pressure_eos;
  MaterialProperty<RankTwoTensor> & _stress_eos_elastic;
  MaterialProperty<RankTwoTensor> & _stress_cpl_elastic;
  MaterialProperty<Real> & _density_corr;
  MaterialProperty<Real> & _mu;

  const Real _x;
  const Real _y;
  MaterialProperty<Real> & _V;
  MaterialProperty<Real> & _pressure_total;


  virtual Real chem_react_1(const Real Y1, const Real T);
  virtual Real chem_react_2(const Real Y1, const Real Y1_rate, const Real Y2, const Real T);
  virtual Real chem_react_3(const Real Y2, const Real Y2_rate, const Real Y3, const Real T);
  virtual Real chem_react_4(const Real Y3, const Real Y3_rate, const Real T);
  usingTensorIndices(i_, j_, k_, l_);
};
