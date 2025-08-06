#include "ADVisHSLUMP.h"

registerMooseObject("catApp", ADVisHSLUMP);

InputParameters
ADVisHSLUMP::validParams()
{
  InputParameters params = ADMatHeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute heating from viscous dissipation");
  params.addRequiredParam<Real>("beta_p", "beta_p");
  params.addRequiredParam<Real>("use_PK2", "use_PK2");
  params.addRequiredCoupledVar("dirac_switch_react", "dirac_switch_react");
  params.addRequiredParam<Real>("thr_activation", "thr_activation");
  return params;
}

ADVisHSLUMP::ADVisHSLUMP(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    //retrieve Ydots
    _cauchy_stress(getMaterialProperty<RankTwoTensor>("cauchy_stress")),
    _beta_p(getParam<Real>("beta_p")),
    _Ep_dot(getMaterialProperty<RankTwoTensor>("Ep_dot")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _F(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _use_PK2(getParam<Real>("use_PK2")),
    _dirac_switch_react(coupledValue("dirac_switch_react")),
    _thr_activation(getParam<Real>("thr_activation"))
{}

ADReal
ADVisHSLUMP::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  ADReal q_plastic;
  
  if (_use_PK2 == 1.){
    Real J = _F[_qp].det();
    RankTwoTensor PK2 = J * _F[_qp].inverse() * _cauchy_stress[_qp] * _F[_qp].inverse().transpose();
    q_plastic = std::max(PK2.doubleContraction(_Ep_dot[_qp]), 0.);
  }else{
    q_plastic = _cauchy_stress[_qp].doubleContraction(_Ep_dot[_qp]);
  }

  ADReal res = _beta_p * (1. / (_rho[_qp] * _cv[_qp])) * q_plastic;

  if(_dirac_switch_react[_qp] > _thr_activation){
    res = _beta_p * (1. / (_rho[_qp] * _cv[_qp])) * q_plastic;
  }else{
    res = 0.;
  }

  return res * _test[_i][_qp];
}