#include "Y2_dot_RDX.h"

registerMooseObject("catApp", Y2_dot_RDX);

InputParameters
Y2_dot_RDX::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription("Y1 species time derivate term of the conservation equation"
                              "this kernel computes only the time derivative contribution"
                              "to the residual");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredCoupledVar("Y1", "Y1");
  return params;
}

Y2_dot_RDX::Y2_dot_RDX(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _r2(getMaterialProperty<Real>("r2")),
    _dr2dT(getMaterialProperty<Real>("dr2dT")),
    _r1(getMaterialProperty<Real>("r1")),
    _dr1dT(getMaterialProperty<Real>("dr1dT")),
    _Y1(coupledValue("Y1")),
    _Y1Id(coupled("Y1")),
    _TempId(coupled("temperature"))
{}

Real
Y2_dot_RDX::computeQpResidual()
{
  Real res;
  res = (_r1[_qp] * _Y1[_qp]) - (_r2[_qp] * _u[_qp]); //Y2_dot
  return - res * _test[_i][_qp];
}

Real 
Y2_dot_RDX::computeQpJacobian()
{
  Real Jac = 0.0;
  Jac = - _r2[_qp];
  return - Jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
Y2_dot_RDX::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  Real OffJac = 0.0;
  if (jvar == _TempId)
    {
    OffJac = (_dr1dT[_qp] * _Y1[_qp]) - (_dr2dT[_qp] * _u[_qp]);
    return - OffJac * _test[_i][_qp] * _phi[_j][_qp];
    }
  else if (jvar == _Y1Id)
    {
    OffJac = _r1[_qp];
    return - OffJac * _test[_i][_qp] * _phi[_j][_qp];
    }
  else
    {
    return 0.0;
    }
}
