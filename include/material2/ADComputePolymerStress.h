#pragma once

#include "ADComputeStressBase.h"
#include "ElasticityTensorTools.h"
#include "Function.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADComputePolymerStress
  : public Material
{
public:
  static InputParameters validParams();

  ADComputePolymerStress(const InputParameters & parameters);
  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const MaterialPropertyName &_bulk_modulus_name;
  const ADMaterialProperty<Real> &_bulk_modulus;
  
  const MaterialPropertyName &_poisson_ratio_name;
  const ADMaterialProperty<Real> &_poisson_ratio;

  const ADMaterialProperty<RankTwoTensor> &_be_bar;
  const ADMaterialProperty<Real> &_J;

  ADMaterialProperty<RankTwoTensor> &_sigma_binder;

private:
  //none
};