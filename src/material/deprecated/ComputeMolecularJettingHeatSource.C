/// Calculates heat generated due to thermal expansion

#include "ComputeMolecularJettingHeatSource.h"

registerMooseObject("TensorMechanicsApp",ComputeMolecularJettingHeatSource);

InputParameters
ComputeMolecularJettingHeatSource::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Moleclar jetting heat source material");
  params.addCoupledVar("c","Phase field damage variable");
  params.addCoupledVar("dirac_switch","dirac delta function to control when the heat source is on/off");
  params.addParam<MaterialPropertyName>("density", "density", "Property name of the density material property");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Property name of the specific_heat material property");
  params.addRequiredParam<Real>("alpha", "Alpha coefficient from curve fit");
  params.addRequiredParam<Real>("up", "Piston velocity");
  params.addRequiredParam<Real>("slope_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("sound_speed","Speed of sound in the material");
  params.addRequiredParam<Real>("gap","Gap dimension in the direction of shock");
  params.addRequiredParam<Real>("zeta","Proportionality constant between gap and up^2");
  params.addRequiredParam<Real>("tau_induction","the time of activation for dirac delta function");
  params.addRequiredParam<Real>("factor","a factor multiplied to control jetting heat");
  return params;
}

ComputeMolecularJettingHeatSource::ComputeMolecularJettingHeatSource(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _c(coupledValue("c")),
    _dirac_switch(coupledValue("dirac_switch")),
    _density(getMaterialProperty<Real>("density")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _alpha(getParam<Real>("alpha")),
    _up(getParam<Real>("up")),
    _s(getParam<Real>("slope_UsUp")),
    _sound_speed(getParam<Real>("sound_speed")),
    _g(getParam<Real>("gap")),
    _zeta(getParam<Real>("zeta")),
    _tau_induction(getParam<Real>("tau_induction")),
    _factor(getParam<Real>("factor")),
    _molecular_jetting_heat(declareProperty<Real>(_base_name + "molecular_jetting_heat")),
    _d_molecular_jetting_heat_dT(declareProperty<Real>(_base_name + "d_molecular_jetting_heat_dT")),
    _d_molecular_jetting_heat_dc(declareProperty<Real>(_base_name + "d_molecular_jetting_heat_dc"))
{
}

void
ComputeMolecularJettingHeatSource::computeQpProperties()
{
  // kbmw = 0.00002806; // kb/mw value in km^2/s^2/K
  Real us = _s * _up + _sound_speed;

  if(_dirac_switch[_qp]>=1. && _dirac_switch[_qp]<=2. ){  
  _molecular_jetting_heat[_qp] = _factor * std::pow(_c[_qp],2.0) * _density[_qp] * _specific_heat[_qp] * _alpha * us *
                                  _up * (1.0 - std::exp(-std::pow(_g *_up / _zeta,1.0) ))/(_tau_induction);

  _d_molecular_jetting_heat_dT[_qp] = 0.0;

  _d_molecular_jetting_heat_dc[_qp] = _factor * 2.0 * _c[_qp] * _density[_qp] * _specific_heat[_qp] * _alpha * us * _up *
                                      (1.0 - std::exp(-std::pow(_g *_up / _zeta,1.0) ))/(_tau_induction);
  }else{
  _molecular_jetting_heat[_qp] = 0.0;

  _d_molecular_jetting_heat_dT[_qp] = 0.0;

  _d_molecular_jetting_heat_dc[_qp] = 0.0;
  }



}
