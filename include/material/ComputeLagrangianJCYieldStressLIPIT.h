
#pragma once

#include "ComputeStressBase.h"
#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "MathUtils.h"
#include "RankFourTensor.h"

class ComputeLagrangianJCYieldStressLIPIT : public Material
{
public:
  static InputParameters validParams();
  ComputeLagrangianJCYieldStressLIPIT(const InputParameters & parameters);

//  void initialSetup() override;

protected:
  const Real _ep_ref;
  MaterialProperty<Real> &_yield_JC;
  const MaterialProperty<Real> &_ep;
  const MaterialProperty<Real> &_ep_old;
  MaterialProperty<Real> &_ep_rate;
  const Real _transtemp;
  const VariableValue &_temperature;
  const Real _T0;
  const Real _k;
  MaterialProperty<Real> &_theta;
  MaterialProperty<Real> &_flow_stress;
  MaterialProperty<Real> &_dH;
  MaterialProperty<Real> &_d2H;
  const Real _A;
  const Real _B;
  const Real _use_temp;
  const Real _a_melt;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  const Real _n_h;
  const Real _use_rate;

  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();
  
};
