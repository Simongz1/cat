#include "ADComputeMixtureYieldStress.h"

registerMooseObject("mistApp", ADComputeMixtureYieldStress);

InputParameters
ADComputeMixtureYieldStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute Johnson-Cook yield stress to be used in lagrangian J2 plasticity");
  params.addRequiredParam<Real>("epsilon_ref", "reference effective plasstic ");

  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("T0", "reference temperature");
  params.addRequiredParam<Real>("k", "temperature exponent");
  params.addRequiredParam<Real>("A", "A");
  params.addRequiredParam<Real>("B", "B");
  params.addRequiredParam<Real>("a_melt", "a_melt");
  params.addRequiredParam<Real>("Tm0", "reference melting temperature");
  params.addRequiredParam<Real>("n_h", "n exponent for hardening");
  //other melting model parameters
  params.addRequiredParam<Real>("melting_model", "melting_model");
  //parameter for mixture
  params.addParam<bool>("use_mixture", "use_mixture");
  params.addCoupledVar("fraction_csv", "fraction_csv");
  params.addParam<Real>("binder_yield", "binder_yield");
  return params;
}

ADComputeMixtureYieldStress::ADComputeMixtureYieldStress(
    const InputParameters & parameters)
  : Material(parameters),
    _ep_ref(getParam<Real>("epsilon_ref")),
    _yield_JC(declareProperty<Real>("yield_JC")),
    _ep(getMaterialProperty<Real>("effective_plastic_strain")), //effective plastic strain from J2 plasticity
    _ep_old(getMaterialPropertyOld<Real>("effective_plastic_strain")), //old effective plastic strain
    _ep_rate(declareProperty<Real>("ep_rate")),
    _temperature(coupledValue("temperature")),
    _T0(getParam<Real>("T0")),
    _k(getParam<Real>("k")),
    _theta(declareProperty<Real>("theta")),
    _flow_stress(declareProperty<Real>("flow_stress")),
    _dH(declareProperty<Real>("dH")),
    _d2H(declareProperty<Real>("d2H")),
    _A(getParam<Real>("A")),
    _B(getParam<Real>("B")),
    _a_melt(getParam<Real>("a_melt")),
    _Tm0(getParam<Real>("Tm0")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _n_h(getParam<Real>("n_h")),
    _melting_temperature(declareProperty<Real>("melting_temperature")),
    _ratedep(declareProperty<Real>("ratedep")),
    _melting_model(getParam<Real>("melting_model")),
    _pressure_total(getADMaterialProperty<Real>("pressure_total")),
    _use_mixture(getParam<bool>("use_mixture")),
    _fraction_csv(coupledValue("fraction_csv")),
    _binder_yield(getParam<Real>("binder_yield"))
{}

void
ADComputeMixtureYieldStress::initQpStatefulProperties()
{}

void
ADComputeMixtureYieldStress::computeQpProperties()
{
  //compute plastic strain rate
  Real ep_rate;
  ep_rate = (_ep[_qp] - _ep_old[_qp]) / _dt;
  _ep_rate[_qp] = ep_rate;

  //melting for HMX
  Real Tm;
  Real J;

  J = _deformation_gradient[_qp].det();

  if (_melting_model == 1){ //lindemann
    Tm = _Tm0 * std::exp(2.0 * _a_melt * (1. - J)) * std::pow((1. / J), 2.0 * (0.7 - _a_melt - 0.33));
  }
  
  if (_melting_model == 0){ //simon
    Tm = _Tm0 * std::pow((1. + (std::abs(_pressure_total[_qp].value()) / 0.9631)), 1. / 2.8855);
  }
  _melting_temperature[_qp] = Tm; //save melting temperature

  Real theta;
  theta = (_temperature[_qp] - _T0) / (Tm - _T0);
  Real captheta = std::min(theta, 0.9); //limit theta to a maximum value of 0.9 to avoid convergence issues

  //compute yield stress value
  Real yield;
  yield = (_A + _B * std::pow(std::max(_ep[_qp], 1e-4), _n_h));
  yield *= (1.0 - std::pow(captheta, _k));

  //include strain rate logarithmic term

  Real ratedep = 1. + (1.086 * std::log(1. + (std::max(ep_rate, 0.) / _ep_ref)));

  _ratedep[_qp] = ratedep;
  _theta[_qp] = captheta;
  _flow_stress[_qp] = _use_mixture ? (_fraction_csv[_qp] * yield * ratedep) + (_binder_yield * (1. - _fraction_csv[_qp])) : 
                    yield * ratedep;
 
  //declare derivatives for hardening convergence
  Real dH;
  Real d2H;

  dH = (_use_mixture ? _fraction_csv[_qp] : 1.) * _B * _n_h * std::pow(std::max(_ep[_qp], 1e-4), _n_h - 1.0);
  d2H = (_use_mixture ? _fraction_csv[_qp] : 1.) * _B * _n_h * (_n_h - 1.0) * std::pow(std::max(_ep[_qp], 1e-4), _n_h - 2.0);

  dH *= (1.0 - std::pow(captheta, _k));
  d2H *= (1.0 - std::pow(captheta, _k));

  _dH[_qp] = dH * ratedep;
  _d2H[_qp] = d2H * ratedep;
}