#include "ADComputeHotspotHeatRate.h"

registerMooseObject("catApp", ADComputeHotspotHeatRate);

InputParameters
ADComputeHotspotHeatRate::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute heat rate for hotspot assignment");
    params.addRequiredCoupledVar("temperature", "temperature variable");
    params.addRequiredCoupledVar("tau_hotspot", "induction time for hotspot initial rate");
    params.addRequiredCoupledVar("target", "target temperature distribution");
    params.addRequiredParam<bool>("use_lump", "whether to compute lumped terms or not");
    params.addRequiredParam<Real>("T_ref", "reference temperature");
    params.addRequiredParam<bool>("direct_T", "apply T directly or no");
    params.addRequiredParam<bool>("use_sin", "ramp HS using a sin");
    params.addRequiredParam<bool>("use_shock", "whether to use shock loading or not");
    params.addRequiredCoupledVar("velocity", "velocity");
    params.addRequiredParam<Real>("element_size", "element_size");
    return params;
}

ADComputeHotspotHeatRate::ADComputeHotspotHeatRate(const InputParameters & parameters)
  : Material(parameters),
    _T(coupledValue("temperature")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _target(coupledValue("target")),
    _tau_hotspot(coupledValue("tau_hotspot")),
    _T_ref(getParam<Real>("T_ref")),
    //declare properties
    _q_hotspot(declareADProperty<Real>("q_hotspot")),
    _use_lump(getParam<bool>("use_lump")),
    _direct_T(getParam<bool>("direct_T")),
    _use_sin(getParam<bool>("use_sin")),
    _use_shock(getParam<bool>("use_shock")),
    _v(coupledValue("velocity")),
    _h(getParam<Real>("element_size"))
{   
}

void
ADComputeHotspotHeatRate::computeQpProperties()
{
    //generic lumped hotspot shape
    ADReal q = (_rho[_qp] * _cv[_qp]) * (_target[_qp] - _T_ref);
    ADReal lump_factor = (1. / (_rho[_qp] * _cv[_qp]));
    
    if(_use_sin){
        ADReal off;
        off = (_t <= _tau_hotspot[_qp] ? 1. : 0.);
        q *= off * (M_PI / (2. * _tau_hotspot[_qp])) * std::sin(M_PI * _t / _tau_hotspot[_qp]);
    }
    
    if(_direct_T){
        q *= 1. / _tau_hotspot[_qp];
    }

    if(_use_lump){
        q *= lump_factor;
    }

    Real window_shock;
    if(_use_shock){
        window_shock = (std::abs(_v[_qp]) > 0.1 && _T[_qp] < _target[_qp] ? 1. : 0.);
        q *= window_shock;
    }

    _q_hotspot[_qp] = q;
}