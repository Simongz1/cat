#include "HMXArrheniusY1.h"

registerMooseObject("beaverApp",HMXArrheniusY1);

InputParameters
HMXArrheniusY1::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Y1 species evolution for Arrhenius model for HMX, values from Tarver 1996");
  params.addRequiredCoupledVar("chemY1", "chemY1");
  params.addRequiredParam<Real>("Z1", "pre-exponential factor for Y1 reaction rate");
  params.addRequiredParam<Real>("E1", "activation energy for reaction 1");
  params.addRequiredParam<Real>("R_const", "gas constant");
  params.addRequiredCoupledVar("temperature", "temperature");
  return params;
}

HMXArrheniusY1::HMXArrheniusY1(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _chemY1(coupledValue("chemY1")),
    _chemY1_dot(getMaterialProperty<Real>("chemY1_dot")),
    _Z1(getParam<Real>("Z1")),
    _E1(getParam<Real>("E1")),
    _R_const(getParam<Real>("R_const")),
    _temperature(coupledValue("temperature"))
{

}

Real
HMXArrheniusY1::computeQpResidual()
{
  //compute the weak form for lambda_dot = f(lambda)

  Real res; //RHS variable for storing the entire form
  res = - _chemY1[_qp] * _Z1 * std::exp(- _E1 / (_R_const * _temperature[_qp]));
  return - res * _test[_i][_qp];
}

Real
HMXArrheniusY1::computeQpJacobian()
{
  //compute the weak form for lambda_dot = f(lambda)

  Real res; //RHS variable for storing the entire form
  res = - _chemY1[_qp] * _Z1 * std::exp(- _E1 / (_R_const * _temperature[_qp]));
  return (res / _chemY1[_qp]) * _phi[_j][_qp] * _test[_i][_qp];
}
