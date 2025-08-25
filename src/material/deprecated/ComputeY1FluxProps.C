#include "ComputeY1FluxProps.h"
#include "RankTwoTensor.h"

registerMooseObject("beaverApp", ComputeY1FluxProps);

InputParameters
ComputeY1FluxProps::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute the required space and time derivatives of density"
                                "velocity and Y1 mass fraction to be used in the Y1 conservation"
                                "equation");
    params.addRequiredCoupledVar("vx", "x component of velocity");
    params.addRequiredCoupledVar("vy", "y component of velocity");
    params.addRequiredCoupledVar("Y1", "Y1 mass fraction");
    return params;
}

ComputeY1FluxProps::ComputeY1FluxProps(const InputParameters & parameters)
  : Material(parameters),
    _density_0(getMaterialProperty<Real>("density_0")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>("mechanica_strain")),
    _epsilon_mech_dot(declareProperty<RankTwoTensor>("epsilon_mech_dot")),
    _vx(coupledValue("vx")),
    _vy(coupledValue("vy")),
    _Y1(coupledValue("Y1")),
    _density_t(declareProperty<Real>("density_t")),
    _density_t_old(getMaterialPropertyOld<Real>("density_t")),
    //the derivatives of the velocities are computed bellow
    _dvx_dx(declareProperty<Real>("dvx_dx")),
    _dvx_dy(declareProperty<Real>("dvx_dy")),
    _dvy_dx(declareProperty<Real>("dvy_dx")),
    _dvy_dy(declareProperty<Real>("dvy_dy")),
    _drho_dx(declareProperty<Real>("drho_dx")),
    _drho_dy(declareProperty<Real>("drho_dy")),
    _drho_dt(declareProperty<Real>("drho_dt")),
    _drho_dt_an(declareProperty<Real>("drho_dt_an")),
    _gradY1(coupledGradient("Y1")),
    _dY1_dx(declareProperty<Real>("dY1_dx")),
    _dY1_dy(declareProperty<Real>("dY1_dy")),
    _gradvx(coupledGradient("vx")),
    _gradvy(coupledGradient("vy"))
{}

void
ComputeY1FluxProps::computeQpProperties()
{
    //define the required derivatives and gradients
    Real dvx_dx = _gradvx[_qp](0);
    Real dvx_dy = _gradvx[_qp](1);
    Real dvy_dx = _gradvy[_qp](0);
    Real dvy_dy = _gradvy[_qp](1);

    //store such gradients in their respective variables
    _dvx_dx[_qp] = dvx_dx;
    _dvx_dy[_qp] = dvx_dy;
    _dvy_dx[_qp] = dvy_dx;
    _dvy_dy[_qp] = dvy_dy;

    _dY1_dx[_qp] = _gradY1[_qp](0);
    _dY1_dy[_qp] = _gradY1[_qp](1);

    //compute density correction and time derivative in terms of mechanical strain
    _density_t[_qp] = _density_0[_qp] / (1.0 + _mechanical_strain[_qp].trace());
    //numerical computation of the time derivative
    _drho_dt[_qp] = (_density_t[_qp] - _density_t_old[_qp]) / _dt;

    //the time derivative can be also analytically derived as
    _epsilon_mech_dot[_qp] = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
    _drho_dt_an[_qp] = _density_0[_qp] * (_epsilon_mech_dot[_qp].trace() / std::pow((1.0 + _mechanical_strain[_qp].trace()), 2.0));

    _drho_dx[_qp] = (1.0 / _vx[_qp]) * _drho_dt[_qp];
    _drho_dy[_qp] = (1.0 / _vy[_qp]) * _drho_dt[_qp];
}