#include "HMXArrheniusY4.h"

registerMooseObject("beaverApp",HMXArrheniusY4);

InputParameters
HMXArrheniusY4::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Y4 species evolution for Arrhenius model for HMX, values from Tarver 1996");
  params.addRequiredCoupledVar("chemY3", "mass fraction of Y4 species in the chemical model");
  params.addRequiredParam<Real>("Z3", "pre-exponential factor for Y4 reaction rate");
  params.addRequiredParam<Real>("E3", "activation energy for reaction 4");
  params.addRequiredParam<Real>("R_const", "gas constant");
  params.addRequiredCoupledVar("temperature", "temperature");
  return params;
}

HMXArrheniusY4::HMXArrheniusY4(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _chemY3(coupledValue("chemY3")),
    _chemY3_dot(getMaterialProperty<Real>("chemY3_dot")),
    _Z3(getParam<Real>("Z3")),
    _E3(getParam<Real>("E3")),
    _R_const(getParam<Real>("R_const")),
    _temperature(coupledValue("temperature"))
{

}

Real
HMXArrheniusY4::computeQpResidual()
{
  //compute the weak form for lambda_dot = f(lambda)

  Real res; //RHS variable for storing the entire form
  res = (_chemY3[_qp] * _chemY3[_qp]) * _Z3 * std::exp(- _E3 / (_R_const * _temperature[_qp]));;
  return - res * _test[_i][_qp];
}

Real
HMXArrheniusY4::computeQpJacobian()
{
  //compute the weak form for lambda_dot = f(lambda)

  Real res; //RHS variable for storing the entire form
  Real jac;
  jac = 0.0;
  return jac * _phi[_j][_qp] * _test[_i][_qp];
}
