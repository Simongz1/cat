#include "ADQYdot.h"

registerMooseObject("beaverApp", ADQYdot);

InputParameters
ADQYdot::validParams()
{
  InputParameters params = ADKernel::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute the heat of reaction contribution to the heat source");
  params.addRequiredParam<Real>("Q1", "Heat of reaction 1");
  params.addRequiredParam<Real>("Q2", "Heat of reaction 2");
  params.addRequiredParam<Real>("Q3", "Heat of reaction 3");
  return params;
}

ADQYdot::ADQYdot(const InputParameters & parameters)
  : ADKernel(parameters),
    _Q1(getParam<Real>("Q1")),
    _Q2(getParam<Real>("Q2")),
    _Q3(getParam<Real>("Q3")),
    //retrieve Ydots as type::ADReal thru getAD
    _Y1dot(getADMaterialProperty<Real>("Y1dot")),
    _Y2dot(getADMaterialProperty<Real>("Y2dot")),
    _Y3dot(getADMaterialProperty<Real>("Y3dot")),
    _Y4dot(getADMaterialProperty<Real>("Y4dot"))
{}

ADReal
ADQYdot::computeQpResidual()
{

  ADReal q_dec_tot = 0.0;
  q_dec_tot += _Q1 * _Y1dot[_qp];
  q_dec_tot += _Q2 * _Y2dot[_qp];
  q_dec_tot += _Q3 * _Y3dot[_qp];
  return q_dec_tot * _test[_i][_qp];

}