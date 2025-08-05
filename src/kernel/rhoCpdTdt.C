#include "rhoCpdTdt.h"

registerMooseObject("beaverApp", rhoCpdTdt);

InputParameters
rhoCpdTdt::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("compute the time derivative of temperature for the heat equation");
  params.addParam<MaterialPropertyName>("Cp", "specific heat material property");
  params.addParam<MaterialPropertyName>("density_t", "name of material property storing time dependent density");
  return params;
}

rhoCpdTdt::rhoCpdTdt(const InputParameters & parameters)
  : DerivativeMaterialInterface<TimeKernel>(parameters),
    _Cp(getMaterialPropertyByName<Real>("Cp")),
    _density_t(getMaterialPropertyByName<Real>("density_t")),
    _u_name(getVar("u", 0) -> name()),
    _dCpdT(getMaterialPropertyDerivative<Real>("dCpdT", _u_name)), //derivative of Cp wrt T
    _drho_dt(getMaterialProperty<Real>("drho_dt"))
{}

Real
rhoCpdTdt::computeQpResidual()
{
    Real res = _density_t[_qp] * _Cp[_qp] * _u_dot[_qp];
    return res * _test[_i][_qp];
}

Real
rhoCpdTdt::computeQpJacobian()
{ 
    //NOTE: for the density derivatives could also use the analytically derive density derivative
    Real Jac = _Cp[_qp] * _drho_dt[_qp] + _density_t[_qp] * _dCpdT[_qp];
    return Jac * _phi[_j][_qp] * _test[_i][_qp];
}