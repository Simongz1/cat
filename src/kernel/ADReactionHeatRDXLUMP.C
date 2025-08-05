#include "ADReactionHeatRDXLUMP.h"

registerMooseObject("catApp", ADReactionHeatRDXLUMP);

InputParameters
ADReactionHeatRDXLUMP::validParams()
{
  InputParameters params = ADMatHeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute the heat of reaction contribution to the heat source");
  params.addRequiredCoupledVar("Y1", "Y1");
  params.addRequiredCoupledVar("Y2", "Y2");
  params.addRequiredCoupledVar("Y3", "Y3");
  params.addRequiredCoupledVar("dirac_switch_react","dirac_switch_react");
  return params;
}

ADReactionHeatRDXLUMP::ADReactionHeatRDXLUMP(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    _Q1(getADMaterialProperty<Real>("Q1")),
    _Q2(getADMaterialProperty<Real>("Q2")),
    _r1(getADMaterialProperty<Real>("r1")),
    _r2(getADMaterialProperty<Real>("r2")),
    _Y1(coupledValue("Y1")),
    _Y2(coupledValue("Y2")),
    _Y3(coupledValue("Y3")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _dirac_switch_react(coupledValue("dirac_switch_react"))
{}

ADReal
ADReactionHeatRDXLUMP::computeQpResidual()
{
  ADReal q_dec_tot = 0.0;
  q_dec_tot -= _Q1[_qp] * (- _r1[_qp] * _Y1[_qp]);
  q_dec_tot += _Q2[_qp] * ((_r2[_qp] * _Y2[_qp]));
  q_dec_tot *= _rho[_qp];

  //control after chemistry heating is done

  ADReal res;
  if (_dirac_switch_react[_qp] > 1){ //here the chemistry heating is done and the continuum model evolves
    res = - (1. / (_rho[_qp] * _cv[_qp])) * q_dec_tot * _test[_i][_qp];
  }else{
    res = 0. * _test[_i][_qp];
  }

  return res;
}