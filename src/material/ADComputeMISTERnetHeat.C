/// Calculates heat generated due to MISTERnet prediction

#include "ADComputeMISTERnetHeat.h"

registerMooseObject("mlApp", ADComputeMISTERnetHeat);

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

  params.addRequiredCoupledVar("v_components", "v_components");
  params.addRequiredCoupledVar("a_components", "a_components");

  params.addRequiredParam<Real>("element_size", "element_size");

  //for surrogate chemistry rate
  params.addRequiredCoupledVar("fraction_csv", "fraction_csv");

  //test: implementation of elemental integral
  params.addParam<bool>("correction_heat", false, "correction_heat: parameter to accont for elemental volume on energy conservation");
  params.addParam<Real>("dirac_tolerance", 0., "dirac_tolerance");
  params.addParam<bool>("use_gating", true, "gating for heat and reaction sources");
  params.addParam<bool>("use_complete_burn", false, "use_complete_burn");
  return params;
}

ADComputeMISTERnetHeat::ADComputeMISTERnetHeat(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _T_ref(getParam<Real>("T_ref")),
    _dirac_switch_shock(adCoupledValue("dirac_switch_shock")),
    _dirac_switch_react(adCoupledValue("dirac_switch_react")),
    _T(adCoupledValue("temperature")),
    _density(getADMaterialProperty<Real>("density")),
    _specific_heat(getADMaterialProperty<Real>("specific_heat")),

    _temperature_mister_shock(getMaterialProperty<Real>("temperature_mister_shock")),
    _temperature_mister_react(getMaterialProperty<Real>("temperature_mister_react")),
    _heatrate_mister_shock(declareADProperty<Real>("heatrate_mister_shock")),
    _heatrate_mister_react(declareADProperty<Real>("heatrate_mister_react")),
    _v_flag(getADMaterialProperty<Real>("v_flag")),

    //Test: formulate a surrogate chemistry evolution source
    _Y1_dot_surrogate(declareADProperty<Real>("Y1_dot_surrogate")),
    _Y2_dot_surrogate(declareADProperty<Real>("Y2_dot_surrogate")),
    _Y3_dot_surrogate(declareADProperty<Real>("Y3_dot_surrogate")),
    _indicator_surrogate(declareProperty<Real>("indicator_surrogate")),

    //for dynamic time update
    _time_react(getADMaterialProperty<Real>("time_react")),
    _time_shock(declareADProperty<Real>("time_shock")),
    _h(getParam<Real>("element_size")),
    _fraction_csv(coupledValue("fraction_csv")),

    _correction_heat(getParam<bool>("correction_heat")),

    //get the shock velocity
    _us(getADMaterialProperty<Real>("us")),
    _use_gating(getParam<bool>("use_gating")),
    _use_complete_burn(getParam<bool>("use_complete_burn"))
{
  const unsigned int n_v = coupledComponents("v_components");
  _v.reserve(n_v);
  for (unsigned int i = 0; i < n_v; ++i)
    _v.push_back(&adCoupledValue("v_components", i));

  const unsigned int n_a = coupledComponents("a_components");
  _a.reserve(n_a);
  for (unsigned int i = 0; i < n_a; ++i)
    _a.push_back(&adCoupledValue("a_components", i));
}

void
ADComputeMISTERnetHeat::computeQpProperties()
{ 
  //obtain velocity and acceleration profiles
  std::vector<ADReal> v_vect(_v.size());
  std::vector<ADReal> a_vect(_a.size());

  for (unsigned int i = 0; i < _v.size(); ++i){
    v_vect[i] = (*_v[i])[_qp];
  }
  for (unsigned int j = 0; j < _a.size(); ++j){
    a_vect[j] = (*_a[j])[_qp];
  }

  //define inline for L2norm
  auto L2norm = [](const std::vector<ADReal> &vect) -> ADReal
  {
    ADReal sum = 0;
    for (int i = 0; i < vect.size(); ++i){
      sum += vect[i] * vect[i];
    }
    return MetaPhysicL::sqrt(sum);
  };

  //define shock and equilibration times
  //TESTING SCALING BY A FACTOR OF 3 IN TIME SHOCK
  ADReal tau_shock = ADReal(_h) / _us[_qp];

  //TESTING: fix the time to not depend on shock time
  //ADReal tau_react = std::max((_time_react[_qp]), tau_shock); 
  ADReal tau_react = std::max(ADReal(_dt), _time_react[_qp]);

  _time_shock[_qp] = tau_shock;

  //define target terms for shock and reaction
  ADReal TargetShock = _density[_qp] * _specific_heat[_qp] * 
                       (_temperature_mister_shock[_qp] - _T_ref);

  ADReal TargetReact = _density[_qp] * _specific_heat[_qp] * 
                       (_temperature_mister_react[_qp] - _temperature_mister_shock[_qp]);

  //define condition for heatrates
  bool condition_shock = _v_flag[_qp] == ADReal(1.) && _dirac_switch_shock[_qp] > 0. && _dirac_switch_shock[_qp] <= _time_shock[_qp] ? true : false;
  bool condition_react = _v_flag[_qp] == ADReal(1.) && _dirac_switch_react[_qp] > 0. && _dirac_switch_react[_qp] <= tau_react ? true : false;

  //assign based on call conditions for shock and reaction
  _heatrate_mister_shock[_qp] = condition_shock ? TargetShock / tau_shock : ADReal(0.);
  _heatrate_mister_react[_qp] = condition_react ? TargetReact / tau_react : ADReal(0.);
  
  bool condition_chem;
  condition_chem = (_temperature_mister_react[_qp] > 1600 ? true : false);

  //construct rates
  //we need to define an indicator to turn on and off the surrogate chemistry source

  if(condition_chem && condition_react){
    _indicator_surrogate[_qp] = 1.;
  }else{
    _indicator_surrogate[_qp] = 0.;
  }

  //define reaction rates for surrogate model
  ADReal Y3_pred_consistent = _use_complete_burn ? ADReal(1.0) : _fraction_csv[_qp];
  ADReal Y3_pred = condition_chem ? Y3_pred_consistent : 0.; //this should go to 1 when complete burn is set
  _Y3_dot_surrogate[_qp] = _indicator_surrogate[_qp] * Y3_pred / tau_react;

  //apply gating
  if (_use_gating){
    _heatrate_mister_shock[_qp] *= CosGate(_time_shock[_qp], _dirac_switch_shock[_qp]);
    _heatrate_mister_react[_qp] *= CosGate(tau_react, _dirac_switch_react[_qp]);
    _Y3_dot_surrogate[_qp] *= CosGate(tau_react, _dirac_switch_react[_qp]);
  }

  //now we have to make Y1 consistent as well
  if (_use_complete_burn){
    //Y1 doesn't start at 1, it starts at _fraction_csv, so we need to account for that in the rate
    _Y1_dot_surrogate[_qp] = - _Y3_dot_surrogate[_qp] * (1. / Y3_pred_consistent) * _fraction_csv[_qp];
  }else{
    _Y1_dot_surrogate[_qp] = - _Y3_dot_surrogate[_qp];
  }
}

ADReal
ADComputeMISTERnetHeat::CosGate(const ADReal tau, const ADReal t){
  //evaluate expression
  ADReal exp = (1. - MetaPhysicL::cos(2 * M_PI * t / tau));
  return exp;
}