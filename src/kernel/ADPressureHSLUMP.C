#include "ADPressureHSLUMP.h"

registerMooseObject("catApp", ADPressureHSLUMP);

InputParameters
ADPressureHSLUMP::validParams()
{
  InputParameters params = ADMatHeatSource::validParams();
  params.addClassDescription("compute elastic thermomechanical coupling heat source term");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredCoupledVar("Yfinal", "reacted mass fraction");
  params.addRequiredParam<Real>("beta_av", "beta_av");
  params.addRequiredCoupledVar("dirac_switch_react", "dirac_switch_react");
  params.addRequiredParam<Real>("thr_activation", "thr_activation");
  return params;
}

ADPressureHSLUMP::ADPressureHSLUMP(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    _T(coupledValue("temperature")),
    _T_var_number(coupled("temperature")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _Yfinal(coupledValue("Yfinal")),
    _pressure_mg(getADMaterialProperty<Real>("pressure_mg")),
    _pressure_JWL(getADMaterialProperty<Real>("pressure_JWL")),
    _beta_av(getParam<Real>("beta_av")),
    _pressure_av(getADMaterialProperty<Real>("pressure_av")),
    _dP_dT(getADMaterialProperty<Real>("dP_dT")),
    _Ee_dot(getMaterialProperty<RankTwoTensor>("Ee_dot")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _dirac_switch_react(coupledValue("dirac_switch_react")),
    _Fe(getMaterialProperty<RankTwoTensor>("Fe")),
    _thr_activation(getParam<Real>("thr_activation"))
{}

ADReal
ADPressureHSLUMP::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  ///full pressure

  ADReal res_pressure;
  RankTwoTensor Ce = _Fe[_qp].transpose() * _Fe[_qp];
  res_pressure = - std::max(_u[_qp] * _dP_dT[_qp] * (Ce.inverse().doubleContraction(_Ee_dot[_qp])), 0.);

  ADReal res_av;
  res_av = _beta_av * _pressure_av[_qp] * _Ee_dot[_qp].trace();

  ADReal res_tot = res_pressure + res_av;

  if(_dirac_switch_react[_qp] > _thr_activation){
    res_tot *= 1.;
  }else{
    res_tot *= 0.;
  }

  return - res_tot * _test[_i][_qp];
}