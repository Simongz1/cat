#include "nablaVY3.h"

registerMooseObject("beaverApp", nablaVY3);

//this kernel computes the divergence term of the conservation
//equation for the chemical reaction species
//as this is a kernel object, the used derivatives are called into the computation
//but explicitly computed within a material object

InputParameters
nablaVY3::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("Y1 species divergence term in conservation equation"
                              "the kernel computes the mass fraction flux divergence"
                              "as $nabla dot (rho Y {vx, vy})$");
  params.addParam<MaterialPropertyName>("density", "variable storing constant denisty");
  params.addCoupledVar("vx", "x component of velocity");
  params.addCoupledVar("vy", "y component of velocity");
  return params;
}

nablaVY3::nablaVY3(const InputParameters & parameters)
  : TimeKernel(parameters),
    _rho(getMaterialPropertyByName<Real>("density")), //time and space dependent density
    _vx(coupledValue("vx")),
    _vy(coupledValue("vy")),
    //coupled gradients
    _grad_vx(coupledGradient("vx")),
    _grad_vy(coupledGradient("vy"))
{}

Real
nablaVY3::computeQpResidual()
{
  //this divergence form is simplified using chain rule tricks
  Real res = 0.0;
  res += _rho[_qp] * (_vx[_qp] * _grad_u[_qp](0) + _u[_qp] * _grad_vx[_qp](0) + _vy[_qp] * _grad_u[_qp](1) + _u[_qp] * _grad_vy[_qp](1));
  return res * _test[_i][_qp];
}
