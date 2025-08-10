#include "ADComputeHotspotHeatRate.h"

registerMooseObject("catApp", ADComputeHotspotHeatRate);

InputParameters
ADComputeHotspotHeatRate::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute heat rate for hotspot assignment");
    params.addRequiredParam<MaterialPropertyName>("target_name", "target temperature distribution name");
    params.addRequiredParam<MaterialPropertyName>("tau_name", "induction time distribution name");
    params.addRequiredParam<bool>("use_lump", "whether to compute lumped terms or not");
    return params;
}

ADComputeHotspotHeatRate::ADComputeHotspotHeatRate(const InputParameters & parameters)
  : Material(parameters),
    
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _target_name(getParam<MaterialPropertyName>("target_name")),
    _target(getADMaterialPropertyByName<Real>(_target_name)),
    _tau_name(getParam<MaterialPropertyName>("tau_name")),
    _tau(getADMaterialProperty<Real>(_tau_name)),

    _use_PK2(getParam<bool>("use_PK2")),
    //declare properties
    _q_hotspot(declareADProperty<Real>("q_hotspot")),
    _use_lump(getParam<bool>("use_lump"))
{   
}

void
ADComputeHotspotHeatRate::computeQpProperties()
{
    //generic lumped hotspot shape
    ADReal q = (1. / _tau[_qp]) * _target[_qp];
    ADReal lump_factor = (1. / (_rho[_qp] * _cv[_qp]));

    if (_use_lump){
        q *= lump_factor;
    }

    _q_hotspot[_qp] = q;
}