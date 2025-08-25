#include "divY1.h"

registerMooseObject("beaverApp", divY1);

//this kernel computes the divergence term of the conservation
//equation for the chemical reaction species
//as this is a kernel object, the used derivatives are called into the computation
//but explicitly computed within a material object

InputParameters
divY1::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Y1 species divergence term in conservation equation"
                              "the kernel computes the mass fraction flux divergence"
                              "as $nabla dot (rho Y {vx, vy})$");
  params.addParam<MaterialPropertyName>("density_t", "variable storing the variable density"
                                        "different from the nominal density rho_0");
  params.addRequiredParam<Real>("vel_dim", "velocity vector dimension");
  params.addCoupledVar("vx", "x component of velocity");
  params.addCoupledVar("vy", "y component of velocity");
  return params;
}

divY1::divY1(const InputParameters & parameters)
  : Kernel(parameters),
    _density_t(getMaterialProperty<Real>("density_t")), //time and space dependent density
    _dvx_dx(getMaterialProperty<Real>("dvx_dt")),
    _dvx_dy(getMaterialProperty<Real>("dvx_dy")),
    _dvy_dx(getMaterialProperty<Real>("dvy_dx")),
    _dvy_dy(getMaterialProperty<Real>("dvy_dy")),
    _drho_dx(getMaterialProperty<Real>("drho_dx")),
    _drho_dy(getMaterialProperty<Real>("drho_dy")),
    _drho_dt_an(getMaterialProperty<Real>("drho_dt_an")),
    _gradY1(getMaterialProperty<Real>("gradY1")),
    _gradvx(getMaterialProperty<Real>("gradvx")),
    _gradvy(getMaterialProperty<Real>("gradvy")),
    _vx(coupledValue("vx")),
    _vy(coupledValue("vy"))
{}

Real
divY1::computeQpResidual()
{
  //initialize divergence
  Real div_x, div_y;
  div_x = (_drho_dx[_qp] * _vx[_qp] * _u[_qp]) + (_density_t[_qp] * _dvx_dx[_qp] * _u[_qp]) + (_density_t[_qp] * _vx[_qp] * _grad_u[_qp](0));
  div_y = (_drho_dy[_qp] * _vy[_qp] * _u[_qp]) + (_density_t[_qp] * _dvy_dy[_qp] * _u[_qp]) + (_density_t[_qp] * _vy[_qp] * _grad_u[_qp](1));
  Real div = div_x + div_y;
  return div * _test[_i][_qp];
}
