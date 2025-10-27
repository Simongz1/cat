
#pragma once

#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "MathUtils.h"
#include "RankFourTensor.h"

class ComputeIntPValRDX : public Material
{
public:
  static InputParameters validParams();

  ComputeIntPValRDX(const InputParameters & parameters);

//  void initialSetup() override;

protected:
  const Real _T_ref;
  const MaterialProperty<Real> &_rho;
  const MaterialProperty<Real> &_Cv;
  const VariableValue &_T;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain_old;
  const Real _C0;
  const Real _C1;

  //variable element length computation
  const MaterialProperty<RankFourTensor> &_elasticity_tensor;
  const Real _Le;
  MaterialProperty<Real> &_pressure_mg;
  MaterialProperty<Real> &_pressure_JWL;
  MaterialProperty<Real> &_pressure_total;
  const VariableValue &_Y_final;
  const Real _A_u;
  const Real _R1_u;
  const Real _B_u;
  const Real _R2_u;
  const Real _omega_u;
  const Real _A_r;
  const Real _R1_r;
  const Real _B_r;
  const Real _R2_r;
  const Real _omega_r;
  MaterialProperty<RankTwoTensor> &_dMG_dT;
  MaterialProperty<RankTwoTensor> &_dJWL_dT;
  const MaterialProperty<Real> &_pressure_total_old;
  MaterialProperty<Real> &_pressure_hist;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> &_f_inv;
  const Real _P0;
  const Real _use_RDX;
  const Real _Gamma;
  const Real _s;
  const Real _A1;
  const Real _R1;
  const Real _A2;
  const Real _R2;
  const Real _omega;
  const VariableValue &_vx;
  const VariableValue &_ax;
  const Real _thr_a;
  const Real _thr_v;
  MaterialProperty<Real> &_v_flag;
  MaterialProperty<Real> &_J;

  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  
};
