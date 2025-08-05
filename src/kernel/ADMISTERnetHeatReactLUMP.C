/// Calculates heat generated due to thermal expansion

#include "ADMISTERnetHeatReactLUMP.h"

registerMooseObject("TensorMechanicsApp", ADMISTERnetHeatReactLUMP);

InputParameters
ADMISTERnetHeatReactLUMP::validParams()
{
  InputParameters params = ADMatHeatSource::validParams();
  params.addClassDescription("Misternet react heat source kernel");
  return params;
}

ADMISTERnetHeatReactLUMP::ADMISTERnetHeatReactLUMP(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    _heatrate_mister_react(getADMaterialProperty<Real>("heatrate_mister_react")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat"))
{}

ADReal
ADMISTERnetHeatReactLUMP::computeQpResidual()
{
  return - (1. / (_rho[_qp] * _cv[_qp])) * _heatrate_mister_react[_qp] * _test[_i][_qp];
}