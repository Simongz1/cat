#include "Y1_dot.h"

registerMooseObject("catApp", Y1_dot);

InputParameters
Y1_dot::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription("create a source term corresponding to the RHS of the conservation equation. Coupled with the tarver model");
  params.addRequiredCoupledVar("temperature", "temperature");
  return params;
}

Y1_dot::Y1_dot(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _r1(getMaterialProperty<Real>("r1")),
    _dr1dT(getMaterialProperty<Real>("dr1dT")),
    _TempId(coupled("temperature"))
{}

Real
Y1_dot::computeQpResidual()
{
  //computes the RHS Ydot term for the chemical species conservation
  //equaton
  //different from dY_dt computation, this directly references a material prop

  Real res;
  res = -1.0 * _r1[_qp] * _u[_qp]; //Y1_dot
  return - res * _test[_i][_qp];
}

Real
Y1_dot::computeQpJacobian()
{
  Real Jac = 0.0;
  //on diagonal jacobian dRYi/dYi
  Jac = -1.0 * _r1[_qp];
  return - Jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
Y1_dot::computeQpOffDiagonalJacobian(unsigned int jvar)
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