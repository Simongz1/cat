  #include "ComputeTempElasticProps.h"

registerMooseObject("TensorMechanicsApp", ComputeTempElasticProps);

InputParameters
ComputeTempElasticProps::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute temperature dependen elastic properties for LIPIT impact model");
  params.addRequiredParam<Real>("epsilon_ref", "reference effective plasstic ");
  params.addRequiredParam<Real>("yield0", "reference yield stress");
  params.addRequiredParam<Real>("yield_brit", "temperature dependent brittle yield stress");
  params.addRequiredParam<Real>("yield_soft", "temperature dependent soft yield stress");
  params.addRequiredParam<Real>("transtemp", "transition temperature");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("T0", "reference temperature");
  params.addRequiredParam<Real>("k", "temperature exponent");
  return params;
}

ComputeTempElasticProps::ComputeTempElasticProps(
    const InputParameters & parameters)
  : Material(parameters),
    _ep_ref(getParam<Real>("epsilon_ref")),
    _yield_JC(declareProperty<Real>("yield_JC")),
    _ep(getMaterialProperty<Real>("effective_plastic_strain")), //effective plastic strain from J2 plasticity
    _ep_old(getMaterialPropertyOld<Real>("effective_plastic_strain")), //old effective plastic strain
    _ep_rate(declareProperty<Real>("ep_rate")),
    _yield0(getParam<Real>("yield0")),
    _yield_brit(getParam<Real>("yield_brit")),
    _yield_soft(getParam<Real>("yield_soft")),
    _transtemp(getParam<Real>("transtemp")),
    _temperature(coupledValue("temperature")),
    _T0(getParam<Real>("T0")),
    _k(getParam<Real>("k")),
    _theta(declareProperty<Real>("theta")),
    _flow_stress(declareProperty<Real>("flow_stress"))
    
{
}

void
ComputeTempElasticProps::computeQpProperties()
{
  //compute plastic strain rate
  Real ep_rate;
  ep_rate = (_ep[_qp] - _ep_old[_qp]) / _dt;
  _ep_rate[_qp] = ep_rate;

  //compute yield stress value
  Real theta;
  theta = (_temperature[_qp] - _T0) / (_transtemp - _T0);
  Real yield;
  yield = _yield0 * std::log(1.0 + (ep_rate / _ep_ref)) * (_yield_brit * std::pow((1.0 - theta), _k) + _yield_soft * std::pow(theta, _k));

  _theta[_qp] = theta;
  _yield_JC[_qp] = yield;
  _flow_stress[_qp] = yield;
}
