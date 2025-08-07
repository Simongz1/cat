/// Calculates heat generated due to MISTERnet prediction

#include "ADComputeMISTERnetHeat.h"

registerMooseObject("TensorMechanicsApp", ADComputeMISTERnetHeat);

InputParameters
ADComputeMISTERnetHeat::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Misternet heat source material, also computes artificial chemistry");
  params.addRequiredParam<Real>("T_ref", "reference temperature for thermal expansion");
  params.addRequiredCoupledVar("dirac_switch_shock","dirac delta function to control when the heat source is on/off");
  params.addRequiredCoupledVar("dirac_switch_react","dirac delta function to control when the heat source is on/off");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addParam<MaterialPropertyName>("density", "density", "Property name of the density material property");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Property name of the specific_heat material property");
  params.addRequiredParam<Real>("factor","a factor multiplied to control jetting heat");
  params.addRequiredParam<Real>("heat_time_shock","time of which heat source is on");
  params.addRequiredParam<Real>("heat_time_react","time of which heat source is on");
  params.addRequiredCoupledVar("vx", "vx");
  params.addRequiredCoupledVar("ax", "ax");
  params.addRequiredParam<Real>("thr_v", "thr_v");
  params.addRequiredParam<Real>("thr_a", "thr_a");
  return params;
}

ADComputeMISTERnetHeat::ADComputeMISTERnetHeat(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _T_ref(getParam<Real>("T_ref")),
    _dirac_switch_shock(coupledValue("dirac_switch_shock")),
    _dirac_switch_react(coupledValue("dirac_switch_react")),
    _T(coupledValue("temperature")),
    _density(getADMaterialProperty<Real>("density")),
    _specific_heat(getADMaterialProperty<Real>("specific_heat")),
    _factor(getParam<Real>("factor")),
    _heat_time_shock(getParam<Real>("heat_time_shock")),
    _heat_time_react(getParam<Real>("heat_time_react")),
    _temperature_mister_shock(getMaterialProperty<Real>("temperature_mister_shock")),
    _temperature_mister_react(getMaterialProperty<Real>("temperature_mister_react")),
    _heatrate_mister_shock(declareADProperty<Real>("heatrate_mister_shock")),
    _heatrate_mister_react(declareADProperty<Real>("heatrate_mister_react")),
    _v_flag(getMaterialProperty<Real>("v_flag")),
    _vx(coupledValue("vx")),
    _ax(coupledValue("ax")),
    _thr_v(getParam<Real>("thr_v")),
    _thr_a(getParam<Real>("thr_a")),

    //Test: formulate a surrogate chemistry evolution source
    _Y1_dot_surrogate(declareADProperty<Real>("Y1_dot_surrogate")),
    _Y2_dot_surrogate(declareADProperty<Real>("Y2_dot_surrogate")),
    _Y3_dot_surrogate(declareADProperty<Real>("Y3_dot_surrogate")),
    _indicator_surrogate(declareProperty<Real>("indicator_surrogate"))

{}

void
ADComputeMISTERnetHeat::computeQpProperties()
{
  

  if(_v_flag[_qp] == 1. && _dirac_switch_shock[_qp] > 0. && _dirac_switch_shock[_qp] < 1.){
    _heatrate_mister_shock[_qp] = std::max((1. / _heat_time_shock) * _density[_qp] * _specific_heat[_qp] * (_temperature_mister_shock[_qp] - _T_ref), 0.);
  }
  else {
    _heatrate_mister_shock[_qp] = 0.0;
  }

  if(_v_flag[_qp] == 1. && _dirac_switch_react[_qp] > 0. && _dirac_switch_react[_qp] < 1.){
    _heatrate_mister_react[_qp] = std::max((1. / _heat_time_react) * _density[_qp] * _specific_heat[_qp] * (_temperature_mister_react[_qp] - _temperature_mister_shock[_qp]), 0.);
  }
  else {
    _heatrate_mister_react[_qp] = 0.0;
  }

  //compute the surrogate rates
  //the surrogate rates depend only on tau_react
  //we generate and assign the predicted values by temperature

  Real Y1_pred;
  Real Y2_pred;
  Real Y3_pred;

  if (_heatrate_mister_react[_qp] > 0.){ //this is the case where deflagration occurs locally
    Y1_pred = 0.;
    Y2_pred = 0.;
    Y3_pred = 1.;
  }else{ //case where q_react < q_shock, meaning quench
    Y1_pred = 1.;
    Y2_pred = 0.;
    Y3_pred = 0.;
  }

  //construct rates
  //we need to define an indicator to turn on and off the surrogate chemistry source

  if(_v_flag[_qp] ==1. && _dirac_switch_react[_qp] > 0. && _dirac_switch_react[_qp] < 1.){ //chemical takeover stage at reacted regions
    _indicator_surrogate[_qp] = 1.;
  }else{
    _indicator_surrogate[_qp] = 0.;
  }

  _Y1_dot_surrogate[_qp] = _indicator_surrogate[_qp] * Y1_pred / _heat_time_react;
  _Y2_dot_surrogate[_qp] = _indicator_surrogate[_qp] * Y2_pred / _heat_time_react;
  _Y3_dot_surrogate[_qp] = _indicator_surrogate[_qp] * Y3_pred / _heat_time_react;
}