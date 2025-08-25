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
  params.addRequiredParam<bool>("direct_T", "assign T predictions directly from a purely numerical source");
  params.addRequiredParam<bool>("dynamic_tau", "use dyanamic reaction time");
  params.addRequiredParam<bool>("temp_crit", "temp_crit");
  params.addRequiredParam<bool>("use_sin", "use_sin");
  params.addRequiredParam<Real>("element_size", "element_size");
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
    _indicator_surrogate(declareProperty<Real>("indicator_surrogate")),
    _direct_T(getParam<bool>("direct_T")),
    //for dynamic time update
    _dynamic_tau(getParam<bool>("dynamic_tau")),
    _time_react(getMaterialProperty<Real>("time_react")),
    _temp_crit(getParam<bool>("temp_crit")),
    _use_sin(getParam<bool>("use_sin")),
    _time_shock(declareProperty<Real>("time_shock")),
    _h(getParam<Real>("element_size"))

{}

void
ADComputeMISTERnetHeat::computeQpProperties()
{
  Real cutoff;
  Real total_tau;
  Real tau_shock;
  if(_dynamic_tau){
    total_tau = std::max(_time_react[_qp], 1e-6);
    cutoff = total_tau;
    tau_shock = _h / std::clamp(std::abs(_vx[_qp]), 0., 10.); //this computes the actual velocity it takes for the shock to cover an element
  }else{
    total_tau = _heat_time_react;
    cutoff = 1.;
    tau_shock = _heat_time_shock;
  }

  _time_shock[_qp] = tau_shock;

  if(_v_flag[_qp] == 1. && _dirac_switch_shock[_qp] > 0. && _dirac_switch_shock[_qp] < 1.){
    if(_direct_T){
      _heatrate_mister_shock[_qp] = _density[_qp] * _specific_heat[_qp] * std::max(_temperature_mister_shock[_qp] - _T_ref, 0.) / _heat_time_shock;
    }else{
      _heatrate_mister_shock[_qp] = std::max((1. / tau_shock) * _density[_qp] * _specific_heat[_qp] * (_temperature_mister_shock[_qp] - _T_ref), 0.);
    }
  }
  else {
    _heatrate_mister_shock[_qp] = 0.0;
  }

  if(_v_flag[_qp] == 1. && _dirac_switch_react[_qp] > 0. && _dirac_switch_react[_qp] < cutoff){
    if(_direct_T){
      _heatrate_mister_react[_qp] = _density[_qp] * _specific_heat[_qp] * std::max(_temperature_mister_react[_qp], 0.) / total_tau;
      if(_use_sin){
        _heatrate_mister_react[_qp] = _density[_qp] * _specific_heat[_qp] * getSinTarget(std::max(_temperature_mister_react[_qp] - _temperature_mister_shock[_qp], 0.), total_tau, std::clamp(total_tau * _dirac_switch_react[_qp], 0., 1.)); 
        //_heatrate_mister_react[_qp] = _density[_qp] * _specific_heat[_qp] * std::max(_temperature_mister_react[_qp], 0.) / total_tau;
        //_heatrate_mister_react[_qp] = _density[_qp] * _specific_heat[_qp] * std::max(_temperature_mister_react[_qp] - _temperature_mister_shock[_qp], 0.) * 
      }
    }else{
      _heatrate_mister_react[_qp] = std::max((1. / total_tau) * _density[_qp] * _specific_heat[_qp] * (_temperature_mister_react[_qp]), 0.);
    }
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

  bool condition;
  if (_temp_crit){
    condition = (_temperature_mister_react[_qp] > 1100. ? true : false);
  }else{
    condition = (_temperature_mister_react[_qp] > _temperature_mister_shock[_qp] ? true : false);
  }

  if (_temperature_mister_react[_qp] > 1100.){ //this is the case where deflagration occurs locally, assuming deflagrated regions occur when T > 1500K
    Y1_pred = 0.;
    Y2_pred = 0.;
    Y3_pred = 1.;
  }else{ //case where q_react < q_shock, meaning quench. We don't do anything in this case
    Y1_pred = 0.;
    Y2_pred = 0.;
    Y3_pred = 0.;
  }

  //construct rates
  //we need to define an indicator to turn on and off the surrogate chemistry source

  if(_v_flag[_qp] == 1. && _dirac_switch_react[_qp] > 0. && _dirac_switch_react[_qp] < cutoff){ //chemical takeover stage at reacted regions
    _indicator_surrogate[_qp] = 1.;
  }else{
    _indicator_surrogate[_qp] = 0.;
  }

  //from conservation of mass: Y1 + Y2 + Y3 = 1

  _Y1_dot_surrogate[_qp] = - _indicator_surrogate[_qp] * Y3_pred / total_tau;
  //_Y2_dot_surrogate[_qp] = _indicator_surrogate[_qp] * Y2_pred / total_tau;
  _Y3_dot_surrogate[_qp] = _indicator_surrogate[_qp] * Y3_pred / total_tau;

  if(_use_sin){
    _Y1_dot_surrogate[_qp] = - _indicator_surrogate[_qp] * Y3_pred / total_tau;
    //_Y1_dot_surrogate[_qp] = - _indicator_surrogate[_qp] * getSinTarget(Y3_pred, total_tau, _dirac_switch_react[_qp]);
    //_Y2_dot_surrogate[_qp] = _indicator_surrogate[_qp] * getSinTarget(Y2_pred, total_tau, _dirac_switch_react[_qp]);
    //_Y3_dot_surrogate[_qp] = _indicator_surrogate[_qp] * getSinTarget(Y3_pred, total_tau, _dirac_switch_react[_qp]);
    _Y3_dot_surrogate[_qp] = _indicator_surrogate[_qp] * Y3_pred / total_tau;
  }
}

Real
ADComputeMISTERnetHeat::getSinTarget(const Real target, const Real induction, const Real time_tracker){
  return target * (M_PI / (2. * induction)) * std::sin(M_PI * time_tracker / induction);
}