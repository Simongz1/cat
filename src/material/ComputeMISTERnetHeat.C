/// Calculates heat generated due to MISTERnet prediction

#include "ComputeMISTERnetHeat.h"

registerMooseObject("TensorMechanicsApp",ComputeMISTERnetHeat);

InputParameters
ComputeMISTERnetHeat::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Misternet heat source material");
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

ComputeMISTERnetHeat::ComputeMISTERnetHeat(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _T_ref(getParam<Real>("T_ref")),
    _dirac_switch_shock(coupledValue("dirac_switch_shock")),
    _dirac_switch_react(coupledValue("dirac_switch_react")),
    _T(coupledValue("temperature")),
    _density(getMaterialProperty<Real>("density")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _factor(getParam<Real>("factor")),
    _heat_time_shock(getParam<Real>("heat_time_shock")),
    _heat_time_react(getParam<Real>("heat_time_react")),
    _temperature_mister_shock(getMaterialProperty<Real>("temperature_mister_shock")),
    _temperature_mister_react(getMaterialProperty<Real>("temperature_mister_react")),
    _heatrate_mister_shock(declareProperty<Real>("heatrate_mister_shock")),
    _heatrate_mister_react(declareProperty<Real>("heatrate_mister_react")),
    _d_heatrate_mister_shock_dT(declareProperty<Real>(_base_name + "_d_heatrate_mister_shock_dT")),
    _d_heatrate_mister_react_dT(declareProperty<Real>(_base_name + "_d_heatrate_mister_react_dT")),
    _v_flag(getMaterialProperty<Real>("v_flag")),
    _vx(coupledValue("vx")),
    _ax(coupledValue("ax")),
    _thr_v(getParam<Real>("thr_v")),
    _thr_a(getParam<Real>("thr_a"))

{
}

//WORKING VERSION
// void
// ComputeMISTERnetHeat::computeQpProperties()
// {
//   //condition to call shock heat source: (abs(vx) > thr)
//   if(_v_flag[_qp]==1. && _dirac_switch_shock[_qp] > 0. && _dirac_switch_shock[_qp] <= 1.){
//     //std::cout << "turned on shock heat" << std::endl;
//     Real heatshock;
//     heatshock = (1.0/_heat_time_shock) * _factor * _density[_qp] * _specific_heat[_qp] * (_temperature_mister_shock[_qp]-_T_ref);
//     if(heatshock <= 0.){
//       _heatrate_mister_shock[_qp]  = 0.;
//       _d_heatrate_mister_shock_dT[_qp] = 0.;
//     }
//     _heatrate_mister_shock[_qp] = heatshock;
//     _d_heatrate_mister_shock_dT[_qp] = _factor * _density[_qp] * _specific_heat[_qp] * (-1.0);
//   }else{
//     _heatrate_mister_shock[_qp] = 0.0;
//     _d_heatrate_mister_shock_dT[_qp] = 0.0;
//   }

//   if(_v_flag[_qp]==1. && _dirac_switch_react[_qp] > 0. && _dirac_switch_react[_qp] <= 1.){
//     //std::cout << "turned on reaction heat" << std::endl;
//     Real heatreact;
//     heatreact = (1.0/_heat_time_react) * _factor * _density[_qp] * _specific_heat[_qp] * (_temperature_mister_react[_qp]-_temperature_mister_shock[_qp]);
//     if (heatreact <= 0.){
//       _heatrate_mister_react[_qp] = 0.;
//       _d_heatrate_mister_react_dT[_qp] = 0.;
//     }
//     _heatrate_mister_react[_qp] = (1.0/_heat_time_react) * _factor * _density[_qp] * _specific_heat[_qp] * (_temperature_mister_react[_qp]-_temperature_mister_shock[_qp]);
//     _d_heatrate_mister_react_dT[_qp] = 0.;
//   }else{
//   _heatrate_mister_react[_qp] = 0.0;
//   _d_heatrate_mister_react_dT[_qp] = 0.0;
//   }
// }

void
ComputeMISTERnetHeat::computeQpProperties()
{
  //the flag already comes turned on from the call conditions
  //originally was the same as react conditions
  if(_v_flag[_qp] == 1. && _dirac_switch_shock[_qp] > 0. ){
    Real heatshock;
    Real qmax = 0.;
    Real MaxDiff = _temperature_mister_shock[_qp] - _T[_qp];
    heatshock = (1. / _heat_time_shock) * _density[_qp] * _specific_heat[_qp] * (_temperature_mister_shock[_qp] - _T[_qp]);
    
    _heatrate_mister_shock[_qp] = heatshock;
    _d_heatrate_mister_shock_dT[_qp] = - (1. / _heat_time_shock) * _density[_qp] * _specific_heat[_qp];

    if (heatshock <= 0.0){
      _heatrate_mister_shock[_qp] = 0.;
      _d_heatrate_mister_shock_dT[_qp] = 0.;
    }
  }
  else {
    _heatrate_mister_shock[_qp] = 0.0;
    _d_heatrate_mister_shock_dT[_qp] = 0.0;
  }

  if(_v_flag[_qp] == 1. && _dirac_switch_react[_qp] > 0. ){
    Real heatreact;
    Real qmax = 0.;
    Real MaxDiff = _temperature_mister_react[_qp] - _T[_qp];
    heatreact = (1. / _heat_time_react) * _density[_qp] * _specific_heat[_qp] * (_temperature_mister_react[_qp] - _T[_qp]);
    
    _heatrate_mister_react[_qp] = heatreact;
    _d_heatrate_mister_react_dT[_qp] = - (1. / _heat_time_react) * _density[_qp] * _specific_heat[_qp];

    if (heatreact <= 0.0){
      _heatrate_mister_react[_qp] = 0.;
      _d_heatrate_mister_react_dT[_qp] = 0.;
    }
  }
  else {
    _heatrate_mister_react[_qp] = 0.0;
    _d_heatrate_mister_react_dT[_qp] = 0.0;
  }
}