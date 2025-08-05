/// Computes EoS stress with artificial viscosity and J2 plasticity with Johnson - Cook model for yield. Radial Return algorithm for returning to the yield surface
/// SGZ 2024

#include "ComputeMieGruneisenChemOnly.h"

registerMooseObject("TensorMechanicsApp", ComputeMieGruneisenChemOnly);

InputParameters
ComputeMieGruneisenChemOnly::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain"
                             "Considers Mie Gruneisen EOS and artificial viscosity damping");
  params.addRequiredParam<Real>("Gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addParam<MaterialPropertyName>("density", "density", "Name of Material Property that provides the density");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Name of Material Property that provides the density");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("slope_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("Le","Maximum element size");

  //chemical reactions params

  params.addRequiredParam<Real>("Z1", "Z1 chemical parameter");
  params.addRequiredParam<Real>("Z2", "Z2 chemical parameter");
  params.addRequiredParam<Real>("Z3", "Z3 chemical parameter");
  params.addRequiredParam<Real>("Z4", "Z4 chemical parameter");

  params.addRequiredParam<Real>("E1", "E1 chemical parameter");
  params.addRequiredParam<Real>("E2", "E2 chemical parameter");
  params.addRequiredParam<Real>("E3", "E3 chemical parameter");

  params.addRequiredParam<Real>("R_const", "gas constant");

  //parameters for JWL EOS

  params.addRequiredParam<Real>("A_g", "A_g parameter for JWL eos");
  params.addRequiredParam<Real>("B_g", "B_g parameter for JWL eos");
  params.addRequiredParam<Real>("C_g", "C_g parameter for JWL eos");
  params.addRequiredParam<Real>("R_1g", "R_1g parameter for JWL eos");
  params.addRequiredParam<Real>("R_2g", "R_2g parameter for JWL eos");
  params.addRequiredParam<Real>("omega", "omega parameter for gas JWL eos");

  //evolution equation parameters

  params.addRequiredParam<Real>("I_ev", "I parameter for evolution equation");
  params.addRequiredParam<Real>("G_ev", "G parameter for evolution equation");
  params.addRequiredParam<Real>("a_ev", "a parameter for evolution equation");
  params.addRequiredParam<Real>("b_ev", "b parameter for evolution equation");
  params.addRequiredParam<Real>("useArrhenius", "chemistry model to use");
  params.addRequiredParam<Real>("x", "x exponent ignition part");
  params.addRequiredParam<Real>("y", "y exponent growth part");
  return params;
}

ComputeMieGruneisenChemOnly::ComputeMieGruneisenChemOnly(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _Gamma(getParam<Real>("Gamma")),
    _density(getMaterialProperty<Real>("density")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _temperature(coupledValue("temperature")),
    _ref_temperature(getParam<Real>("reference_temperature")),
    _s(getParam<Real>("slope_UsUp")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("Le")),
    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _ss_prop(declareProperty<Real>(_base_name + "ss_prop")),
    //chemical reaction variables
    _chemY1(declareProperty<Real>(_base_name + "chemY1")),
    _chemY2(declareProperty<Real>(_base_name + "chemY2")),
    _chemY3(declareProperty<Real>(_base_name + "chemY3")),
    _chemY4(declareProperty<Real>(_base_name + "chemY4")),
    _chemY1_dot(declareProperty<Real>(_base_name + "chemY1_dot")),
    _chemY2_dot(declareProperty<Real>(_base_name + "chemY2_dot")),
    _chemY3_dot(declareProperty<Real>(_base_name + "chemY3_dot")),
    _chemY4_dot(declareProperty<Real>(_base_name + "chemY4_dot")),
    _Z1(getParam<Real>("Z1")),
    _Z2(getParam<Real>("Z2")),
    _Z3(getParam<Real>("Z3")),
    _Z4(getParam<Real>("Z4")),
    _E1(getParam<Real>("E1")),
    _E2(getParam<Real>("E2")),
    _E3(getParam<Real>("E3")),
    _R_const(getParam<Real>("R_const")),
    _chemY1_old(getMaterialPropertyOld<Real>(_base_name + "chemY1")),
    _chemY2_old(getMaterialPropertyOld<Real>(_base_name + "chemY2")),
    _chemY3_old(getMaterialPropertyOld<Real>(_base_name + "chemY3")),
    _chemY4_old(getMaterialPropertyOld<Real>(_base_name + "chemY4")),
    //gas phase parameters
    _A_g(getParam<Real>("A_g")),
    _B_g(getParam<Real>("B_g")),
    _C_g(getParam<Real>("C_g")),
    _R_1g(getParam<Real>("R_1g")),
    _R_2g(getParam<Real>("R_2g")),
    _omega(getParam<Real>("omega")),
    _p_g(declareProperty<Real>("p_g")),
    _p_s(declareProperty<Real>("p_s")),
    // chemical decomposition evolution equation parameters
    _lambda_dot(declareProperty<Real>("lambda_dot")),
    _lambda(declareProperty<Real>("lambda")),
    _lambda_old(getMaterialPropertyOld<Real>("lambda")),
    _I_ev(getParam<Real>("I_ev")),
    _G_ev(getParam<Real>("G_ev")),
    _a_ev(getParam<Real>("a_ev")),
    _b_ev(getParam<Real>("b_ev")),
    _useArrhenius(getParam<Real>("useArrhenius")),
    _pressure_eos(declareProperty<Real>("pressure_eos")),
    _stress_eos_elastic(declareProperty<RankTwoTensor>("stress_eos_elastic")),
    _stress_cpl_elastic(declareProperty<RankTwoTensor>("stress_cpl_elastic")),
    _density_corr(declareProperty<Real>("density_corr")),
    _mu(declareProperty<Real>("mu")),
    _x(getParam<Real>("x")),
    _y(getParam<Real>("y")),
    _V(declareProperty<Real>("V")),
    _pressure_total(declareProperty<Real>("pressure_tota"))

    
{
}

void
ComputeMieGruneisenChemOnly::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();

  //initialize the chemical composition

  _chemY1[_qp] = 1.0;
  _chemY2[_qp] = 0.0;
  _chemY3[_qp] = 0.0;
  _chemY4[_qp] = 0.0;
  
  //_chemY1_dot[_qp] = 0.0;
  //_chemY2_dot[_qp] = 0.0;s
  //_chemY3_dot[_qp] = 0.0;
  //_chemY4_dot[_qp] = 0.0;

  _lambda[_qp] = 0.0;
 
}

void
ComputeMieGruneisenChemOnly::computeQpStress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  //  Calculate "elastic" stresses due to equation of state
  // Calculate pressure from Mie Gruneisen (Menon, 2014), (Zhang, 2011)
  Real K0, delta, eta, temperature, p_s, ss_prop;
  RankTwoTensor stress_eos, stress, stress_cpl;
  K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  ss_prop = std::sqrt(K0 / _density[_qp]); //store sound speed as a function of actual bulk modulus from elasticity tensor
  _ss_prop[_qp] = ss_prop;
  _bulk_modulus[_qp] = K0;
  delta = _mechanical_strain[_qp].trace();
  eta = - delta;
  temperature = _temperature[_qp];

  Real V = 1.0 + _mechanical_strain[_qp].trace(); //relative volume change
  _V[_qp] = V; //store variable
  _density_corr[_qp] = (1.0 / V) * _density[_qp]; //density correction
  _mu[_qp] = - delta / (1.0 + delta);

  //eos for the unreacte phase (Y1): Mie-Gruneisen
  p_s = _Gamma * _density[_qp] * _specific_heat[_qp] * (_ref_temperature - temperature) * (1 / V) + (((K0 * eta) / std::pow((1.0 - (_s * eta)), 2.0)) * ((0.5 * _Gamma) * ((1.0 / V) - 1.0) - 1.0));
  _p_s[_qp] = - p_s; //write into the saved material property variable
  
  //eos for the reacted phases (1 - Y1 = Y2 + Y3 + Y4): Jones-Wilkins-Lee C_form

  _p_g[_qp] = _A_g * std::exp(- _R_1g * (1.0 - eta)) + _B_g * std::exp(- _R_2g * (1.0 - eta)) + _C_g * std::pow((1.0 - eta), - (1.0 + _omega));

  //combined phase eos depending on chemistry model
  Real peos;
  if (_useArrhenius == 1) {
    peos = ((1.0 - _chemY1[_qp]) * _p_g[_qp]) + (_chemY1[_qp] * _p_s[_qp]);
  }
  else {
    peos = ((1.0 - _lambda[_qp]) * _p_s[_qp]) + (_lambda[_qp] * _p_g[_qp]);
  }

  // Calculate bulk-viscosity stress term
  Real trD, jacob, q_bv;
  trD = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  jacob = 1.0 + _mechanical_strain[_qp].trace();
  q_bv = 0.0;
  if (jacob < 1.0) {
    q_bv = ( _C0 * _density[_qp] * trD * std::abs(trD) * std::pow(_Le,2.0) / std::pow(jacob,2.0) ) + ( _C1 * _density[_qp] * ss_prop * trD * _Le / jacob ); //edited sound speed CHECK
  }
  
  _pressure_total[_qp] = peos + q_bv; 

  //evolution equation for chemistry single step reaction

  _lambda_dot[_qp] = (_I_ev * std::pow((1.0 - _lambda[_qp]), _x) * std::pow(_mu[_qp], _a_ev)) + (_G_ev * (1.0 - _lambda[_qp]) * std::pow(_lambda[_qp], _y) * std::pow(_pressure_total[_qp], _b_ev));
  _lambda[_qp] = _lambda_old[_qp] + (_lambda_dot[_qp] * _dt);

  //update variables
  _pressure_eos[_qp] = peos;
  stress_eos = peos * I2;
  _stress_eos_elastic[_qp] = stress_eos;
  stress_cpl = _elasticity_tensor[_qp] * (_elastic_strain_old[_qp] + _strain_increment[_qp]) - K0 * delta * I2;
  _stress_cpl_elastic[_qp] = stress_cpl;
  stress = stress_eos + stress_cpl;
  _stress[_qp] = stress;
  _stress[_qp] += q_bv * I2;

  //chemical reaction evolution functions
  if (_useArrhenius == 1) {
    _chemY1_dot[_qp] = chem_react_1(_chemY1[_qp], _temperature[_qp]);
    _chemY2_dot[_qp] = chem_react_2(_chemY1[_qp], _chemY1_dot[_qp], _chemY2[_qp], _temperature[_qp]);
    _chemY3_dot[_qp] = chem_react_3(_chemY2[_qp], _chemY2_dot[_qp], _chemY3[_qp], _temperature[_qp]);
    _chemY4_dot[_qp] = chem_react_4(_chemY3[_qp], _chemY3_dot[_qp], _temperature[_qp]);

    //formulate the explicit integration to obtain each species

    _chemY1[_qp] = _chemY1_old[_qp] + (_chemY1_dot[_qp] * _dt);
    _chemY2[_qp] = _chemY2_old[_qp] + (_chemY2_dot[_qp] * _dt);
    _chemY3[_qp] = _chemY3_old[_qp] + (_chemY3_dot[_qp] * _dt);
    _chemY4[_qp] = _chemY4_old[_qp] + (_chemY4_dot[_qp] * _dt);
  }
}

Real
ComputeMieGruneisenChemOnly::chem_react_1(const Real Y1, const Real T) {
  Real Y1_rate;
  Y1_rate = - Y1 * _Z1 * std::exp(- (_E1) / (_R_const * T));
  return Y1_rate;
}

Real
ComputeMieGruneisenChemOnly::chem_react_2(const Real Y1, const Real Y1_rate, const Real Y2, const Real T) {
  Real Y2_rate;
  Y2_rate = (Y1 * _Z1 * std::exp(- (_E1) / (_R_const * T))) - (Y2 * _Z2 * std::exp(- (_E2) / (_R_const * T)));
  return Y2_rate;
}

Real 
ComputeMieGruneisenChemOnly::chem_react_3(const Real Y2, const Real Y2_rate, const Real Y3, const Real T) {
  Real Y3_rate;
  Y3_rate = (Y2 * _Z2 * std::exp(- (_E2) / (_R_const * T))) - ((Y3 * Y3) * _Z3 * std::exp(- (_E3) / (_R_const * T)));
  return Y3_rate;
}

Real 
ComputeMieGruneisenChemOnly::chem_react_4(const Real Y3, const Real Y3_rate, const Real T) {
  Real Y4_rate;
  Y4_rate = (Y3 * Y3) * _Z3 * std::exp(- (_E3) / (_R_const * T));
  return Y4_rate;
}