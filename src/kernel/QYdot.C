#include "QYdot.h"

registerMooseObject("catApp", QYdot);

InputParameters
QYdot::validParams()
{
  InputParameters params = HeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute the heat of reaction contribution to the heat source");
  params.addRequiredParam<Real>("Q1", "Heat of reaction 1");
  params.addRequiredParam<Real>("Q2", "Heat of reaction 2");
  params.addRequiredParam<Real>("Q3", "Heat of reaction 3");
  params.addRequiredCoupledVar("Y1", "Y1");
  params.addRequiredCoupledVar("Y2", "Y2");
  params.addRequiredCoupledVar("Y3", "Y3");
  params.addRequiredCoupledVar("Y4", "Y4");
  return params;
}

QYdot::QYdot(const InputParameters & parameters)
  : HeatSource(parameters),
    _Q1(getParam<Real>("Q1")),
    _Q2(getParam<Real>("Q2")),
    _Q3(getParam<Real>("Q3")),
    //retrieve Ydots
    _Y1dot(getMaterialProperty<Real>("Y1dot")),
    _Y2dot(getMaterialProperty<Real>("Y2dot")),
    _Y3dot(getMaterialProperty<Real>("Y3dot")),
    _Y4dot(getMaterialProperty<Real>("Y4dot")),
    //retrieve rate derivatives wrt temperature
    _r1(getMaterialProperty<Real>("r1")),
    _r2(getMaterialProperty<Real>("r2")),
    _r3(getMaterialProperty<Real>("r3")),
    _dr1dT(getMaterialProperty<Real>("dr1dT")),
    _dr2dT(getMaterialProperty<Real>("dr2dT")),
    _dr3dT(getMaterialProperty<Real>("dr3dT")),
    //coupled variables
    _Y1(coupledValue("Y1")),
    _grad_Y1(coupledGradient("Y1")),
    _Y1Id(coupled("Y1")),
    _Y2(coupledValue("Y2")),
    _grad_Y2(coupledGradient("Y2")),
    _Y2Id(coupled("Y2")),
    _Y3(coupledValue("Y3")),
    _grad_Y3(coupledGradient("Y3")),
    _Y3Id(coupled("Y3")),
    _rho(getMaterialProperty<Real>("density"))
{}

Real
QYdot::computeQpResidual()
{

  Real q_dec_tot = 0.0;
  q_dec_tot += _Q1 * (- _r1[_qp] * _Y1[_qp]);
  q_dec_tot += _Q2 * ((_r1[_qp] * _Y1[_qp]) - (_r2[_qp] * _Y2[_qp]));
  q_dec_tot += _Q3 * ((_r2[_qp] * _Y2[_qp]) - (_r3[_qp] * std::pow(_Y3[_qp], 2.0)));
  return - q_dec_tot * _test[_i][_qp];

}

Real
QYdot::computeQpJacobian()
{
  Real Jac = 0.0;
  Jac += _Q1 * (- _dr1dT[_qp] * _Y1[_qp]);
  Jac += _Q2 * ((_dr1dT[_qp] * _Y1[_qp]) - (_dr2dT[_qp] * _Y2[_qp]));
  Jac += _Q3 * ((_dr2dT[_qp] * _Y2[_qp]) - (_dr3dT[_qp] * std::pow(_Y3[_qp], 2.0)));
  return - Jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
QYdot::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  Real OffJac;
  if (jvar == _Y1Id)
    {
    OffJac = _Q1 * (- _r1[_qp]) + _Q2 * (_r1[_qp]);
    }
  else if (jvar == _Y2Id)
    {
    OffJac = _Q2 * (- _r2[_qp]) + _Q3 * (_r2[_qp]);
    }
  else if (jvar == _Y3Id)
    {
    OffJac = - 2.0 * _r3[_qp] * _Y3[_qp];
    }
  else
    {
    OffJac = 0.0;
    }
  return - OffJac * _test[_i][_qp] * _phi[_j][_qp];
}