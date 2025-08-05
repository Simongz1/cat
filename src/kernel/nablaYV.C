#include "nablaYV.h"

registerMooseObject("beaverApp", nablaYV);

//this kernel computes the divergence term of the conservation
//equation for the chemical reaction species
//as this is a kernel object, the used derivatives are called into the computation
//but explicitly computed within a material object

InputParameters
nablaYV::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("kernel to compute the second term of the evolution equation for the chemical species conservation");
  params.addParam<MaterialPropertyName>("density", "variable storing constant denisty");
  params.addCoupledVar("vx", "x component of velocity");
  params.addCoupledVar("vy", "y component of velocity");
  return params;
}

nablaYV::nablaYV(const InputParameters & parameters)
  : TimeKernel(parameters),
    _rho(getMaterialProperty<Real>("density")), //time and space dependent density
    _vx(coupledValue("vx")),
    _vy(coupledValue("vy")),
    //coupled gradients
    _grad_vx(coupledGradient("vx")),
    _grad_vy(coupledGradient("vy"))
{}

Real
nablaYV::computeQpResidual()
{
  //this divergence form is simplified using chain rule tricks
  Real res = 0.0;
  res = (_rho[_qp]) * (_vx[_qp] * _grad_u[_qp](0) + _u[_qp] * _grad_vx[_qp](0) + _vy[_qp] * _grad_u[_qp](1) + _u[_qp] * _grad_vy[_qp](1));
  return res * _test[_i][_qp];
}

Real
nablaYV::computeQpJacobian()
{
  Real Jac = 0.0;
  return Jac;
}
