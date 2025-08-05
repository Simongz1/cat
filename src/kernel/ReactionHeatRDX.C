#include "ReactionHeatRDX.h"

registerMooseObject("catApp", ReactionHeatRDX);

InputParameters
ReactionHeatRDX::validParams()
{
  InputParameters params = HeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute the heat of reaction contribution to the heat source");
  params.addRequiredCoupledVar("Y1", "Y1");
  params.addRequiredCoupledVar("Y2", "Y2");
  params.addRequiredCoupledVar("Y3", "Y3");
  return params;
}

ReactionHeatRDX::ReactionHeatRDX(const InputParameters & parameters)
  : HeatSource(parameters),
    _Q1(getMaterialProperty<Real>("Q1")),
    _Q2(getMaterialProperty<Real>("Q2")),
    //retrieve Ydots
    _Y1dot(getMaterialProperty<Real>("Y1dot")),
    _Y2dot(getMaterialProperty<Real>("Y2dot")),
    _Y3dot(getMaterialProperty<Real>("Y3dot")),
    //retrieve rate derivatives wrt temperature
    _r1(getMaterialProperty<Real>("r1")),
    _r2(getMaterialProperty<Real>("r2")),
    _dr1dT(getMaterialProperty<Real>("dr1dT")),
    _dr2dT(getMaterialProperty<Real>("dr2dT")),
    //coupled variables
    _Y1(coupledValue("Y1")),
    _Y1Id(coupled("Y1")),
    _Y2(coupledValue("Y2")),
    _Y2Id(coupled("Y2")),
    _Y3(coupledValue("Y3")),
    _Y3Id(coupled("Y3")),
    _rho(getMaterialProperty<Real>("density"))
{}

Real
ReactionHeatRDX::computeQpResidual()
{

  Real q_dec_tot = 0.0;
  q_dec_tot -= _Q1[_qp] * (- _r1[_qp] * _Y1[_qp]);
  q_dec_tot += _Q2[_qp] * ((_r2[_qp] * _Y2[_qp]));
  return - q_dec_tot * _test[_i][_qp];

}

Real
ReactionHeatRDX::computeQpJacobian()
{
  Real Jac = 0.0;
  Jac += _Q1[_qp] * (- _dr1dT[_qp] * _Y1[_qp]);
  Jac += _Q2[_qp] * ((_dr1dT[_qp] * _Y1[_qp]) - (_dr2dT[_qp] * _Y2[_qp]));
  return - Jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
ReactionHeatRDX::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  return - 0.0 * _test[_i][_qp] * _phi[_j][_qp];
}