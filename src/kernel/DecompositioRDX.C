#include "DecompositionRDX.h"

registerMooseObject("beaverApp", DecompositionRDX);

InputParameters
DecompositionRDX::validParams()
{
  InputParameters params = HeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute the heat of reaction contribution to the heat source");
  params.addRequiredCoupledVar("Y1", "Y1");
  params.addRequiredCoupledVar("Y2", "Y2");
  params.addRequiredCoupledVar("Y3", "Y3");
  return params;
}

DecompositionRDX::DecompositionRDX(const InputParameters & parameters)
  : HeatSource(parameters),
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
    _rho(getMaterialProperty<Real>("density")),
    _Q1(getMaterialProperty<Real>("Q1")),
    _Q2(getMaterialProperty<Real>("Q2"))
{}

Real
DecompositionRDX::computeQpResidual()
{

  Real q_dec_tot = 0.0;
  q_dec_tot += _Q1[_qp] * (- _r1[_qp] * _Y1[_qp]);
  q_dec_tot += _Q2[_qp] * (_r2[_qp] * _Y2[_qp]);
  return - _rho[_qp] * q_dec_tot * _test[_i][_qp];

}

Real
DecompositionRDX::computeQpJacobian()
{
  Real Jac = 0.0;
  Jac += _Q1[_qp] * (- _dr1dT[_qp] * _Y1[_qp]);
  Jac += _Q2[_qp] * (_dr2dT[_qp] * _Y2[_qp]);
  return - _rho[_qp] * Jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
DecompositionRDX::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  Real OffJac;
  if (jvar == _Y1Id)
    OffJac = - _Q1[_qp] * (- _r1[_qp]);
  else if (jvar == _Y2Id)
    OffJac = 0.0;
  else if (jvar == _Y3Id)
    OffJac = _Q2[_qp] * _r2[_qp];
  else
    OffJac = 0.0;
  return OffJac * _test[_i][_qp] * _phi[_j][_qp];
}