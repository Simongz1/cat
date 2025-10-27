/// Calculates heat generated due to molecular jetting

#pragma once

#include "HeatSource.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// This kernel calculates the heat source term corresponding to molecular jetting

class MISTERnetHeatReact : public HeatSource
{
public:
  static InputParameters validParams();

  MISTERnetHeatReact(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  // virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  std::string _base_name;


  const MaterialProperty<Real> & _heatrate_mister_shock;
  const MaterialProperty<Real> & _heatrate_mister_react;

  const MaterialProperty<Real> & _d_heatrate_mister_shock_dT;
  const MaterialProperty<Real> & _d_heatrate_mister_react_dT;
};
