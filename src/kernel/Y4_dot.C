#include "Y4_dot.h"

registerMooseObject("catApp", Y4_dot);

InputParameters
Y4_dot::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription("Y1 species time derivate term of the conservation equation"
                              "this kernel computes only the time derivative contribution"
                              "to the residual");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredCoupledVar("Y3", "Y3");
  return params;
}

Y4_dot::Y4_dot(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _r3(getMaterialProperty<Real>("r3")),
    _dr3dT(getMaterialProperty<Real>("dr3dT")),
    _Y3(coupledValue("Y3")),
    _Y3Id(coupled("Y3")),
    _TempId(coupled("temperature"))
{}

Real
Y4_dot::computeQpResidual()
{
  Real res;
  res = (_r3[_qp] * std::pow(_Y3[_qp], 2.0)); //Y1_dot
  return - res * _test[_i][_qp];
}

Real 
Y4_dot::computeQpJacobian()
{
  Real Jac = 0.0;
  Jac = 0.0;
  return - Jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
Y4_dot::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  Real OffJac = 0.0;
  if (jvar == _TempId)
    {
    OffJac = (_dr3dT[_qp] * std::pow(_Y3[_qp], 2.0));
    return - OffJac * _test[_i][_qp] * _phi[_j][_qp];
    }
  else if (jvar == _Y3Id)
    {
    OffJac = 2.0 * _r3[_qp] * _Y3[_qp];
    return - OffJac * _test[_i][_qp] * _phi[_j][_qp];
    }
  else
    {
    return 0.0;
    }
}
