/// Calculates heat generated due to thermal expansion

#include "ADMISTERnetHeatShock.h"

registerMooseObject("TensorMechanicsApp", ADMISTERnetHeatShock);

InputParameters
ADMISTERnetHeatShock::validParams()
{
  InputParameters params = ADMatHeatSource::validParams();
  params.addClassDescription("Misternet react heat source kernel");
  return params;
}

ADMISTERnetHeatShock::ADMISTERnetHeatShock(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    _heatrate_mister_shock(getADMaterialProperty<Real>("heatrate_mister_shock"))
{}

ADReal
ADMISTERnetHeatShock::computeQpResidual()
{
  return - _heatrate_mister_shock[_qp] * _test[_i][_qp];
}