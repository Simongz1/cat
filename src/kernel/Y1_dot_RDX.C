#include "Y1_dot_RDX.h"

registerMooseObject("catApp", Y1_dot_RDX);

InputParameters
Y1_dot_RDX::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription("compute tarver model reaction rate for Y1");
  params.addRequiredCoupledVar("temperature", "temperature");
  return params;
}

Y1_dot_RDX::Y1_dot_RDX(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _r1(getMaterialProperty<Real>("r1")),
    _dr1dT(getMaterialProperty<Real>("dr1dT")),
    _TempId(coupled("temperature"))
{}

Real
Y1_dot_RDX::computeQpResidual()
{
  Real res;
  res = -1.0 * _r1[_qp] * _u[_qp]; //- r1 * Y1
  return - res * _test[_i][_qp];
}

Real
Y1_dot_RDX::computeQpJacobian()
{
  Real Jac = 0.0;
  //on diagonal jacobian
  Jac = -1.0 * _r1[_qp];
  return - Jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
Y1_dot_RDX::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  Real OffJac = 0.0;
  if (jvar == _TempId)
    {
    OffJac = -1.0 * _dr1dT[_qp] * _u[_qp];
    return - OffJac * _test[_i][_qp] * _phi[_j][_qp];
    }
  else
    {
    return 0.0;
    }
}