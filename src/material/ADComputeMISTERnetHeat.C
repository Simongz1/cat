/// Calculates heat generated due to MISTERnet prediction

#include "ADComputeMISTERnetHeat.h"

registerMooseObject("mistApp", ADComputeMISTERnetHeat);

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

  params.addRequiredCoupledVar("v_vect", "v_vect");
  params.addRequiredCoupledVar("a_vect", "a_vect");

  params.addRequiredParam<bool>("temp_crit", "temp_crit");
  params.addRequiredParam<Real>("element_size", "element_size");

  //for surrogate chemistry rate
  params.addRequiredCoupledVar("fraction_csv", "fraction_csv");

  //test: implementation of elemental integral
  params.addParam<bool>("correction_heat", false, "correction_heat: parameter to accont for elemental volume on energy conservation");
  params.addParam<Real>("dirac_tolerance", 0., "dirac_tolerance");
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

    _temperature_mister_shock(getMaterialProperty<Real>("temperature_mister_shock")),
    _temperature_mister_react(getMaterialProperty<Real>("temperature_mister_react")),
    _heatrate_mister_shock(declareADProperty<Real>("heatrate_mister_shock")),
    _heatrate_mister_react(declareADProperty<Real>("heatrate_mister_react")),
    _v_flag(getMaterialProperty<Real>("v_flag")),
    _v_vect(coupledVectorValue("v_vect")),
    _a_vect(coupledVectorValue("a_vect")),


    //Test: formulate a surrogate chemistry evolution source
    _Y1_dot_surrogate(declareADProperty<Real>("Y1_dot_surrogate")),
    _Y2_dot_surrogate(declareADProperty<Real>("Y2_dot_surrogate")),
    _Y3_dot_surrogate(declareADProperty<Real>("Y3_dot_surrogate")),
    _indicator_surrogate(declareProperty<Real>("indicator_surrogate")),
    //for dynamic time update
    _time_react(getMaterialProperty<Real>("time_react")),
    _temp_crit(getParam<bool>("temp_crit")),
    _time_shock(declareProperty<Real>("time_shock")),
    _h(getParam<Real>("element_size")),
    _fraction_csv(coupledValue("fraction_csv")),

    _dirac_tolerance(getParam<Real>("dirac_tolerance")),
    _correction_heat(getParam<bool>("correction_heat"))

{}

void
ADComputeMISTERnetHeat::computeQpProperties()
{

  Real tau_react = std::max(_time_react[_qp], _dt);
  Real tau_shock = _h / std::clamp(_v_vect[_qp].norm(), 1., 10.); //this computes the actual velocity it takes for the shock to cover an element

  _time_shock[_qp] = tau_shock;

  if(_v_flag[_qp] == 1. && _dirac_switch_shock[_qp] > 0. && _dirac_switch_shock[_qp] < 1.){
    _heatrate_mister_shock[_qp] = (_density[_qp] * _specific_heat[_qp] * 
                                  (_temperature_mister_shock[_qp] - _T_ref));
  }
  else {
    _heatrate_mister_shock[_qp] = 0.0;
  }

  //test
  ADReal energy_react = _density[_qp] * _specific_heat[_qp] * (_temperature_mister_react[_qp] - _temperature_mister_shock[_qp]);

  const unsigned int npoints = _qrule->n_points();
  
  ADReal correction = _correction_heat ? npoints : 1.;
  ADReal q_dot = correction * energy_react;

  if(_v_flag[_qp] == 1. && _dirac_switch_react[_qp] > 0. && _dirac_switch_react[_qp] < 1. + _dirac_tolerance){
    _heatrate_mister_react[_qp] = q_dot;
  }
  else {
    _heatrate_mister_react[_qp] = 0.0;
  }

  bool condition;
  if (_temp_crit){
    condition = (_temperature_mister_react[_qp] > 1000 ? true : false);
  }else{
    condition = (_temperature_mister_react[_qp] > _temperature_mister_shock[_qp] ? true : false);
  }

  //construct rates
  //we need to define an indicator to turn on and off the surrogate chemistry source

  if(condition && _v_flag[_qp] == 1. && _dirac_switch_react[_qp] > 0. && _dirac_switch_react[_qp] < 1. + _dirac_tolerance){ //reaction window
    _indicator_surrogate[_qp] = 1.;
  }else{
    _indicator_surrogate[_qp] = 0.;
  }

  const Real Y3_pred = (_heatrate_mister_react[_qp] > 0.) ? _fraction_csv[_qp] : 0.;

  //from conservation of mass: Y1 + Y2 + Y3 = 1

  _Y1_dot_surrogate[_qp] = - _indicator_surrogate[_qp] * Y3_pred;
  _Y3_dot_surrogate[_qp] = - _Y1_dot_surrogate[_qp];
}

Real
ADComputeMISTERnetHeat::getSinTarget(const Real target, const Real induction, const Real time_tracker){
  const Real phase = std::clamp(time_tracker / induction, 0., 1.);
  return (M_PI * target / (2. * induction)) * std::sin(M_PI * phase);
}