/// Calculates heat generated due to thermal expansion

#include "ADMISTERnetHeatReact.h"

registerMooseObject("TensorMechanicsApp", ADMISTERnetHeatReact);

InputParameters
ADMISTERnetHeatReact::validParams()
{
  InputParameters params = ADMatHeatSource::validParams();
  params.addClassDescription("Misternet react heat source kernel");
  return params;
}

ADMISTERnetHeatReact::ADMISTERnetHeatReact(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    _heatrate_mister_react(getADMaterialProperty<Real>("heatrate_mister_react"))
{}

ADReal
ADMISTERnetHeatReact::computeQpResidual()
{
  return - _heatrate_mister_react[_qp] * _test[_i][_qp];
}