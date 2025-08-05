#include "ComputeLagrangianJCYieldStressLIPIT.h"

registerMooseObject("TensorMechanicsApp", ComputeLagrangianJCYieldStressLIPIT);

InputParameters
ComputeLagrangianJCYieldStressLIPIT::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute Johnson-Cook yield stress to be used in lagrangian J2 plasticity");
  params.addRequiredParam<Real>("epsilon_ref", "reference effective plasstic ");
  params.addRequiredParam<Real>("transtemp", "transition temperature");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("T0", "reference temperature");
  params.addRequiredParam<Real>("k", "temperature exponent");
  params.addRequiredParam<Real>("A", "A");
  params.addRequiredParam<Real>("B", "B");
  params.addRequiredParam<Real>("use_temp", "use_temp");
  params.addRequiredParam<Real>("a_melt", "a_melt");
  params.addRequiredParam<Real>("n_h", "n exponent for hardening");
  params.addRequiredParam<Real>("use_rate", "use strain");
  return params;
}

ComputeLagrangianJCYieldStressLIPIT::ComputeLagrangianJCYieldStressLIPIT(
    const InputParameters & parameters)
  : Material(parameters),
    _ep_ref(getParam<Real>("epsilon_ref")),
    _yield_JC(declareProperty<Real>("yield_JC")),
    _ep(getMaterialProperty<Real>("effective_plastic_strain")), //effective plastic strain from J2 plasticity
    _ep_old(getMaterialPropertyOld<Real>("effective_plastic_strain")), //old effective plastic strain
    _ep_rate(declareProperty<Real>("ep_rate")),
    _transtemp(getParam<Real>("transtemp")),
    _temperature(coupledValue("temperature")),
    _T0(getParam<Real>("T0")),
    _k(getParam<Real>("k")),
    _theta(declareProperty<Real>("theta")),
    _flow_stress(declareProperty<Real>("flow_stress")),
    _dH(declareProperty<Real>("dH")),
    _d2H(declareProperty<Real>("d2H")),
    _A(getParam<Real>("A")),
    _B(getParam<Real>("B")),
    _use_temp(getParam<Real>("use_temp")),
    _a_melt(getParam<Real>("a_melt")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _n_h(getParam<Real>("n_h")),
    _use_rate(getParam<Real>("use_rate"))
    
{
}

void
ComputeLagrangianJCYieldStressLIPIT::initQpStatefulProperties()
{
	//_flow_stress[_qp] = _yield_brit; //this initializes the yield stress value to compute first time step
}

void
ComputeLagrangianJCYieldStressLIPIT::computeQpProperties()
{
  //compute plastic strain rate
  Real ep_rate;
  ep_rate = (_ep[_qp] - _ep_old[_qp]) / _dt;
  _ep_rate[_qp] = ep_rate;

  //melting for HMX
  Real Tm;
  Real J;

  J = _deformation_gradient[_qp].det();
  //Tm = _Tm0 * std::exp(2.0 * _a_melt * J) * std::pow(J, - 2.0 * (0.7 - _a_melt - 0.33));

  Real theta;
  theta = (_temperature[_qp] - _T0) / (_transtemp - _T0);
  Real captheta = std::min(theta, 0.95); //limit theta to a maximum value of 0.9 to avoid convergence issues

  //compute yield stress value
  Real yield;
  yield = (_A + _B * std::pow(std::max(_ep[_qp], 1e-4), _n_h)); //* (1.0 + std::log(1.0 + (_ep[_qp] / _ep_ref))) * (1.0 - std::pow(theta, _k));
  if (_use_temp == 1){
    yield *= (1.0 - std::pow(captheta, _k));
  }

  //include strain rate logarithmic term

  Real ratedep;
  ratedep = (1.0 + std::log(1.0 + (std::max(0., ep_rate) / _ep_ref)));

  _theta[_qp] = captheta;

  _flow_stress[_qp] = yield;

  if(_use_rate == 1.0){
    _flow_stress[_qp] *= ratedep;
  }

  //declare derivatives for hardening convergence
  Real dH;
  Real d2H;

  dH = _B * _n_h * std::pow(std::max(_ep[_qp], 1e-4), _n_h - 1.0);
  d2H = _B * _n_h * (_n_h - 1.0) * std::pow(std::max(_ep[_qp], 1e-4), _n_h - 2.0);

  if(_use_temp == 1){
    dH *= (1.0 - std::pow(captheta, _k));
    d2H *= (1.0 - std::pow(captheta, _k));
  }
  //write into the material property for derivatives retrieved
  if(_use_rate == 1.0){
    _dH[_qp] = dH * ratedep;
    _d2H[_qp] = d2H * ratedep;
  }
}
