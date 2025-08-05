/// Calculates heat generated due to thermal expansion

#include "MISTERnetHeatReact.h"

registerMooseObject("TensorMechanicsApp",MISTERnetHeatReact);

InputParameters
MISTERnetHeatReact::validParams()
{
  InputParameters params = HeatSource::validParams();
  params.addClassDescription("Misternet react heat source kernel");
  return params;
}

MISTERnetHeatReact::MISTERnetHeatReact(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _heatrate_mister_shock(getMaterialProperty<Real>("heatrate_mister_shock")),
    _heatrate_mister_react(getMaterialProperty<Real>("heatrate_mister_react")),
    _d_heatrate_mister_shock_dT(getMaterialProperty<Real>(_base_name + "_d_heatrate_mister_shock_dT")),
    _d_heatrate_mister_react_dT(getMaterialProperty<Real>(_base_name + "_d_heatrate_mister_react_dT"))
{
}

Real
MISTERnetHeatReact::computeQpResidual()
{
  return - _heatrate_mister_react[_qp] * _test[_i][_qp];
}

Real
MISTERnetHeatReact::computeQpJacobian()
{
  return - _d_heatrate_mister_react_dT[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}


