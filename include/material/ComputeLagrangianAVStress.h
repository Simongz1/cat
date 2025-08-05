
#pragma once

#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "MathUtils.h"
#include "RankFourTensor.h"

class ComputeLagrangianAVStress : public Material
{
public:
  static InputParameters validParams();
  ComputeLagrangianAVStress(const InputParameters & parameters);
  
protected:
  virtual void computeQpProperties() override;
  const MaterialProperty<Real> &_rho;
  const Real _C0;
  const Real _C1;
  const Real _Le;
  MaterialProperty<RankTwoTensor> &_pressure_av;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient_old;
  const MaterialProperty<RankFourTensor> &_elasticity_tensor;
};
