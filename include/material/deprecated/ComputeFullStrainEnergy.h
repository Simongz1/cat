#pragma once

#include "Material.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "PiecewiseBilinear.h"
#include "RankFourTensor.h"

class ComputeFullStrainEnergy : public Material
{
public:
  static InputParameters validParams();
  ComputeFullStrainEnergy(const InputParameters & parameters);


protected:
  const MaterialProperty<Real> &_alpha;
  MaterialProperty<Real> &_W_tot;
  MaterialProperty<Real> &_W_el;
  MaterialProperty<Real> &_W_th;
  MaterialProperty<Real> &_W_bk;

  MaterialProperty<Real> &_dW_bk_dIv1;
  MaterialProperty<Real> &_dW_bk_dIv2;
  MaterialProperty<Real> &_dW_bk_dIv3;

  MaterialProperty<Real> &_Iv1;
  MaterialProperty<Real> &_Iv2;
  MaterialProperty<Real> &_Iv3;

  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient_old;
  const MaterialProperty<RankFourTensor> &_elasticity_tensor;

  const Real _T0;
  const VariableValue &_temperature;
  const Real _nu;
  const Real _beta;
  MaterialProperty<RankTwoTensor> &_cauchy_bk;
  MaterialProperty<RankTwoTensor> &_lagrangian_strain;
  MaterialProperty<RankTwoTensor> &_lagrangian_strain_rate;

  virtual void computeQpProperties() override;
};

