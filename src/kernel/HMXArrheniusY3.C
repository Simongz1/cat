#include "HMXArrheniusY3.h"

registerMooseObject("beaverApp",HMXArrheniusY3);

InputParameters
HMXArrheniusY3::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Y3 species evolution for Arrhenius model for HMX, values from Tarver 1996");
  params.addRequiredCoupledVar("chemY3", "mass fraction of Y3 species in the chemical model");
  params.addRequiredParam<Real>("Z3", "pre-exponential factor for Y3 reaction rate");
  params.addRequiredParam<Real>("E3", "activation energy for reaction 3");
  params.addRequiredParam<Real>("R_const", "gas constant");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("Z2", "Z2");
  params.addRequiredParam<Real>("E2", "E2");
  params.addRequiredCoupledVar("chemY2", "chemY2");
  return params;
}

HMXArrheniusY3::HMXArrheniusY3(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _chemY3(coupledValue("chemY3")),
    _chemY3_dot(getMaterialProperty<Real>("chemY3_dot")),
    _Z3(getParam<Real>("Z3")),
    _E3(getParam<Real>("E3")),
    _R_const(getParam<Real>("R_const")),
    _temperature(coupledValue("temperature")),
    _Z2(getParam<Real>("Z2")),
    _E2(getParam<Real>("E2")),
    _chemY2(coupledValue("chemY2"))
{

}

Real
HMXArrheniusY3::computeQpResidual()
{
  //compute the weak form for lambda_dot = f(lambda)

  Real res; //RHS variable for storing the entire form
  res = _chemY2[_qp] * _Z2 * std::exp(- _E2 / (_R_const * _temperature[_qp]));
  res += - (_chemY3[_qp] * _chemY3[_qp]) * _Z3 * std::exp(- _E3 / (_R_const * _temperature[_qp]));
  return - res * _test[_i][_qp];
}

Real
HMXArrheniusY3::computeQpJacobian()
{
  //compute the weak form for lambda_dot = f(lambda)

  Real res; //RHS variable for storing the entire form
  Real jac;
  jac = - 2.0 * _chemY3[_qp] * _Z3 * std::exp(- _E3 / (_R_const * _temperature[_qp]));
  return jac * _phi[_j][_qp] * _test[_i][_qp];
}
