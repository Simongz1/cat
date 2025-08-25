
#pragma once

#include "ComputeStressBase.h"
#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "MathUtils.h"
#include "RankFourTensor.h"

class ComputeTempElasticProps : public Material
{
public:
  static InputParameters validParams();

  ComputeTempElasticProps(const InputParameters & parameters);

//  void initialSetup() override;

protected:
  const Real _ep_ref;
  MaterialProperty<Real> &_yield_JC;
  const MaterialProperty<Real> &_ep;
  const MaterialProperty<Real> &_ep_old;
  MaterialProperty<Real> &_ep_rate;
  const Real _yield0;
  const Real _yield_brit;
  const Real _yield_soft;
  const Real _transtemp;
  const VariableValue &_temperature;
  const Real _T0;
  const Real _k;
  MaterialProperty<Real> &_theta;
  MaterialProperty<Real> &_flow_stress;

  virtual void computeQpProperties();
};
