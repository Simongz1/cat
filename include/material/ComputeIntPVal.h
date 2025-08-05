
#pragma once

#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "MathUtils.h"
#include "RankFourTensor.h"

class ComputeIntPVal : public Material
{
public:
  static InputParameters validParams();

  ComputeIntPVal(const InputParameters & parameters);

//  void initialSetup() override;

protected:
  const Real _Gamma;
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
  MaterialProperty<Real> &_Je;
  const Real _s;
  const VariableValue &_Y_final;
  const Real _A1;
  const Real _R1;
  const Real _A2;
  const Real _R2;
  const Real _omega;
  MaterialProperty<Real> &_dMG_dT;
  MaterialProperty<Real> &_dJWL_dT;
  MaterialProperty<Real> &_pressureHS;
  MaterialProperty<Real> &_pressureHS2;
  MaterialProperty<Real> &_pressure_flag;
  const MaterialProperty<Real> &_pressure_total_old;
  const Real _flag_threshold;
  MaterialProperty<Real> &_pressure_hist;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient_old;
  MaterialProperty<RankTwoTensor> &_lagrangian_strain;

  virtual void computeQpProperties() override;
};
