#include "ADNablaRhoVY.h"

registerMooseObject("beaverApp", ADNablaRhoVY);

//this kernel computes the divergence term of the conservation
//equation for the chemical reaction species
//as this is a kernel object, the used derivatives are called into the computation
//but explicitly computed within a material object

InputParameters
ADNablaRhoVY::validParams()
{
  InputParameters params = ADTimeKernel::validParams();
  params.addClassDescription("implement the general contribution to the residual for the term nabla rho V Y_i");
  params.addParam<MaterialPropertyName>("density", "variable storing constant denisty");
  params.addCoupledVar("vx", "x component of velocity");
  params.addCoupledVar("vy", "y component of velocity");
  return params;
}

ADNablaRhoVY::ADNablaRhoVY(const InputParameters & parameters)
  : ADTimeKernel(parameters),
    _rho(getADMaterialPropertyByName<Real>("density")), //time and space dependent density
    _vx(adCoupledValue("vx")),
    _vy(adCoupledValue("vy")),
    //coupled gradients
    _grad_vx(adCoupledGradient("vx")),
    _grad_vy(adCoupledGradient("vy"))
{}

ADReal
ADNablaRhoVY::computeQpResidual()
{
  //this divergence form is simplified using chain rule tricks
  ADReal res = 0.0;
  res += _rho[_qp] * (_vx[_qp] * _grad_u[_qp](0) + _u[_qp] * _grad_vx[_qp](0) + _vy[_qp] * _grad_u[_qp](1) + _u[_qp] * _grad_vy[_qp](1));
  return res * _test[_i][_qp];
}
