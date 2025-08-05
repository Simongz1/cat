#include "Y3_dot_RDX.h"

registerMooseObject("catApp", Y3_dot_RDX);

InputParameters
Y3_dot_RDX::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription("Y1 species time derivate term of the conservation equation"
                              "this kernel computes only the time derivative contribution"
                              "to the residual");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredCoupledVar("Y2", "Y2");
  return params;
}

Y3_dot_RDX::Y3_dot_RDX(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _r2(getMaterialProperty<Real>("r2")),
    _dr2dT(getMaterialProperty<Real>("dr2dT")),
    _Y2(coupledValue("Y2")),
    _Y2Id(coupled("Y2")),
    _TempId(coupled("temperature"))
{}

Real
Y3_dot_RDX::computeQpResidual()
{
  Real res;
  res = (_r2[_qp] * _Y2[_qp]); //Y1_dot
  return - res * _test[_i][_qp];
}

Real 
Y3_dot_RDX::computeQpJacobian()
{
  Real Jac = 0.0;
  Jac = 0.0;
  return - Jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
Y3_dot_RDX::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  Real OffJac = 0.0;
  if (jvar == _TempId)
    {
    OffJac = (_dr2dT[_qp] * _Y2[_qp]);
    return - OffJac * _test[_i][_qp] * _phi[_j][_qp];
    }
  else if (jvar == _Y2Id)
    {
    OffJac = _r2[_qp];
    return - OffJac * _test[_i][_qp] * _phi[_j][_qp];
    }
  else
    {
    return 0.0;
    }
}
