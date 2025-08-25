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
    params.addRequiredCoupledVar("dirac_switch_react", "dirac_switch_react");
    params.addRequiredParam<bool>("use_lump", "whether to compute lumped terms or not");
    params.addRequiredParam<Real>("T_ref", "reference temperature");
    params.addRequiredParam<bool>("direct_T", "apply T directly or no");
    params.addRequiredParam<bool>("use_current", "use current T difference");
    params.addParam<Real>("k_proportional", "proportionality constant for temperature");
    params.addRequiredParam<bool>("use_sin", "ramp HS using a sin");
    params.addRequiredParam<Real>("pi", "number pi");
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
    _use_current(getParam<bool>("use_current")),
    _k_proportional(getParam<Real>("k_proportional")),
    _dirac_switch_react(coupledValue("dirac_switch_react")),
    _use_sin(getParam<bool>("use_sin")),
    _pi(getParam<Real>("pi"))
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
    
    if(_use_lump){
        q *= lump_factor;
    }

    if (_t > 0.05 && _t < 0.1) { // around 0.31 ns, first qp
    mooseInfo("t=", _t,
                " tau=", _tau_hotspot[_qp],
                " target=", _target[_qp],
                " T_ref=", _T_ref,
                " dT=", _target[_qp] - _T_ref,
                " q_spec(K/s)=", _q_hotspot[_qp]);
    }

    _q_hotspot[_qp] = q;
}