#pragma once

#include "ADMaterial.h"
//#include "ElasticityTensorTools.h"
//#include "ADSingleVariableReturnMappingSolution.h"
//#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADComputeDeformationGradient
  : public DerivativeMaterialInterface<ADMaterial>//,
    //public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADComputeDeformationGradient(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const ADMaterialProperty<RankTwoTensor> &_strain_increment;
  const ADMaterialProperty<RankTwoTensor> &_rotation_increment;
  ADMaterialProperty<RankTwoTensor> &_F;
  const MaterialProperty<RankTwoTensor> &_F_old;
  ADMaterialProperty<RankTwoTensor> &_C;
  ADMaterialProperty<RankTwoTensor> &_epsilon;
};