/// Calculates heat generated due to thermal expansion

#include "ADMISTERnetHeatShockLUMP.h"

registerMooseObject("TensorMechanicsApp", ADMISTERnetHeatShockLUMP);

InputParameters
ADMISTERnetHeatShockLUMP::validParams()
{
  InputParameters params = ADMatHeatSource::validParams();
  params.addClassDescription("Misternet react heat source kernel");
  return params;
}

ADMISTERnetHeatShockLUMP::ADMISTERnetHeatShockLUMP(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    _heatrate_mister_shock(getADMaterialProperty<Real>("heatrate_mister_shock")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat"))
{}

ADReal
ADMISTERnetHeatShockLUMP::computeQpResidual()
{
  return - (1. / (_rho[_qp] * _cv[_qp])) * _heatrate_mister_shock[_qp] * _test[_i][_qp];
}