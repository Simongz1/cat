#include "ADVisHS.h"

registerMooseObject("catApp", ADVisHS);

InputParameters
ADVisHS::validParams()
{
  InputParameters params = ADMatHeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute heating from viscous dissipation");
  params.addRequiredParam<Real>("beta_p", "beta_p");
  return params;
}

ADVisHS::ADVisHS(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    //retrieve Ydots
    _cauchy_stress(getMaterialProperty<RankTwoTensor>("cauchy_stress")),
    _beta_p(getParam<Real>("beta_p")),
    _Ep_dot(getMaterialProperty<RankTwoTensor>("Ep_dot"))
{}

ADReal
ADVisHS::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  ADReal q_plastic = _cauchy_stress[_qp].doubleContraction(_Ep_dot[_qp]);
  return _beta_p * q_plastic * _test[_i][_qp];
}