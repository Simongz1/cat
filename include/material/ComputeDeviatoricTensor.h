#pragma once

#include "Material.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "PiecewiseBilinear.h"
#include "RankFourTensor.h"
#include "RankTwoScalarTools.h"

class ComputeDeviatoricTensor : public Material
{
public:
  static InputParameters validParams();
  ComputeDeviatoricTensor(const InputParameters & parameters);


protected:
  const   MaterialPropertyName _property_name;
  const MaterialProperty<RankTwoTensor> & _tensorFull;
  MaterialProperty<RankTwoTensor> & _tensorDev;

  virtual void computeQpProperties() override;
};

