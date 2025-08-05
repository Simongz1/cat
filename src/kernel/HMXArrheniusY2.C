#include "HMXArrheniusY2.h"

registerMooseObject("beaverApp",HMXArrheniusY2);

InputParameters
HMXArrheniusY2::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Y2 species evolution for Arrhenius model for HMX, values from Tarver 1996");
  params.addRequiredCoupledVar("chemY2", "mass fraction of Y2 species in the chemical model");
  params.addRequiredParam<Real>("Z2", "pre-exponential factor for Y2 reaction rate");
  params.addRequiredParam<Real>("E2", "activation energy for reaction 2");
  params.addRequiredParam<Real>("R_const", "gas constant");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("Z1", "Z1");
  params.addRequiredParam<Real>("E1", "E1");
  params.addRequiredCoupledVar("chemY1", "chemY1");
  return params;
}

HMXArrheniusY2::HMXArrheniusY2(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _chemY2(coupledValue("chemY2")),
    _chemY2_dot(getMaterialProperty<Real>("chemY2_dot")),
    _Z2(getParam<Real>("Z2")),
    _E2(getParam<Real>("E2")),
    _R_const(getParam<Real>("R_const")),
    _temperature(coupledValue("temperature")),
    _Z1(getParam<Real>("Z1")),
    _E1(getParam<Real>("E1")),
    _chemY1(coupledValue("chemY1"))
{

}

Real
HMXArrheniusY2::computeQpResidual()
{
  //compute the weak form for lambda_dot = f(lambda)

  Real res; //RHS variable for storing the entire form
  res = _chemY1[_qp] * _Z1 * std::exp(- _E1 / (_R_const * _temperature[_qp]));
  res += - _chemY2[_qp] * _Z2 * std::exp(- _E2 / (_R_const * _temperature[_qp]));
  return - res * _test[_i][_qp];
}

Real
HMXArrheniusY2::computeQpJacobian()
{
  //compute the weak form for lambda_dot = f(lambda)

  Real res; //RHS variable for storing the entire form
  Real jac;
  jac = - _Z2 * std::exp(- _E2 / (_R_const * _temperature[_qp]));
  return jac * _phi[_j][_qp] * _test[_i][_qp];
}
