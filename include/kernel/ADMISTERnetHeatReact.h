/// Calculates heat generated due to molecular jetting

#pragma once

#include "ADMatHeatSource.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// This kernel calculates the heat source term corresponding to molecular jetting

class ADMISTERnetHeatReact : public ADMatHeatSource
{
public:
  static InputParameters validParams();

  ADMISTERnetHeatReact(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();
  const ADMaterialProperty<Real> & _heatrate_mister_react;
};
