#include "EosHeatSourceNew.h"

registerMooseObject("catApp",EosHeatSourceNew);

InputParameters
EosHeatSourceNew::validParams()
{
  //InputParameters params = validParams<HeatSource>();
  InputParameters params = HeatSource::validParams();
  params.addClassDescription("Thermal expansion heat source kernel"
                             "generic kernel for finite strain"
                             "for Mie Gruneisen equation of state"
                             "(Menon, 2014) (Zhang, 2011)");
  // params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addRequiredParam<Real>("gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  // params.addRequiredParam<Real>("s_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  // params.addRequiredParam<Real>("K0", "reference bulk modulus");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  params.addRequiredParam<Real>("h_max","Maximum element size");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  // params.addRequiredParam<Real>("ss", "sound speed");
  params.addRequiredParam<Real>("beta_av", "artificial viscosity beta parameter");
  params.addRequiredParam<Real>("beta_p", "plastic flow beta parameter");
  return params;
}

EosHeatSourceNew::EosHeatSourceNew(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _gamma(getParam<Real>("gamma")), // Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS
    // _K0(getParam<Real>("K0")),
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density_name")),
    _deformation_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient
    _deformation_gradient_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient, previous timestep
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")), //elasticity tensor
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _elastic_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _total_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _total_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "total_strain")),
    _h_max(getParam<Real>("h_max")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    // _ss(getParam<Real>("ss")),
    _beta_av(getParam<Real>("beta_av")),
    _beta_p(getParam<Real>("beta_p")),
    _plastic_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_strain")),
    _ss_prop(getMaterialProperty<Real>(_base_name + "ss_prop")) //comes from ComputeMieGruneisenPlastic...
{

}

Real
EosHeatSourceNew::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  Real heat_source;
  Real thermal_expansion_coeff, q_eos; // thermal expansion coefficient depends on Gruneisen parameter, bulk modulus and sound speed
  RankTwoTensor epsilon_T, epsilon_e, epsilon_p, epsilon_m, epsilon_T_dot, epsilon_e_dot, epsilon_p_dot, epsilon_m_dot;

  // compute total, elastic and plastic infinitesimal strains with time derivatives
  epsilon_T = _total_strain[_qp];
  epsilon_e = _elastic_strain[_qp];
  epsilon_p = _plastic_strain[_qp]; // from e_T = e_e + e_p
  epsilon_m = _mechanical_strain[_qp];

  epsilon_T_dot = (_total_strain[_qp] - _total_strain_old[_qp]) / _dt;
  epsilon_e_dot = (_elastic_strain[_qp] - _elastic_strain_old[_qp]) / _dt;
  epsilon_p_dot = (_plastic_strain[_qp] - _plastic_strain_old[_qp]) / _dt; // by conventional derivative additive property
  epsilon_m_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  
  Real epsilon_T_trace_dot, epsilon_e_trace_dot, epsilon_p_trace_dot, epsilon_m_trace_dot;
  epsilon_T_trace_dot = (_total_strain[_qp].trace() - _total_strain_old[_qp].trace()) / _dt;
  epsilon_e_trace_dot =	(_elastic_strain[_qp].trace() - _elastic_strain_old[_qp].trace()) / _dt;
  epsilon_p_trace_dot = (_plastic_strain[_qp].trace() - _plastic_strain_old[_qp].trace()) / _dt;
  epsilon_m_trace_dot = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;

  // compute artificial viscosity stress and heat source terms
  Real q_tot;
  q_tot = 0.0;

  Real stress_av, q_av, J;
  // stress_av = _density[_qp] * _C0 * std::pow(_h_max, 2) * (std::pow(epsilon_T_dot.trace(), 2) - _C1 * _ss_prop[_qp] * epsilon_T_dot.trace());
  J = 1 + _mechanical_strain[_qp].trace();
  stress_av = ( _C0 * _density[_qp] * epsilon_m_trace_dot * std::abs(epsilon_m_trace_dot) * std::pow(_h_max,2.0) / std::pow(J,2.0) ) + ( _C1 * _density[_qp] * _ss_prop[_qp] * epsilon_m_trace_dot * _h_max / J );
  q_av = _beta_av * stress_av * epsilon_m_dot.trace();
  
  if (J < 1.0) {
      q_tot += std::abs(q_av);  // to ensure that only when in compression the av term contributes to heat generation
  }
  
  // compute plastic flow stress and heat source
  // Real q_p;
  RankTwoTensor stress_e, q_p;
  stress_e = _elasticity_tensor[_qp] * (_mechanical_strain[_qp] - _plastic_strain[_qp]); //elastic only stress
  q_p = _beta_p * stress_e * epsilon_p_dot;
  q_tot += std::abs(q_p.doubleContraction(I2));

  // compute thermomechanical coupled term
  Real temp, e0, alpha;
  temp = _u[_qp]; // temperature variable

  q_eos = - temp * _gamma * _density[_qp] * _specific_heat[_qp] * epsilon_m_dot.trace();
  q_tot += std::abs(q_eos);

  return - q_tot * _test[_i][_qp];
}

Real
EosHeatSourceNew::computeQpJacobian()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  Real heat_source;
  Real thermal_expansion_coeff, q_eos; // thermal expansion coefficient depends on Gruneisen parameter, bulk modulus and sound speed
  RankTwoTensor epsilon_T, epsilon_e, epsilon_p, epsilon_m, epsilon_T_dot, epsilon_e_dot, epsilon_p_dot, epsilon_m_dot;

  // compute total, elastic and plastic infinitesimal strains with time derivatives
  epsilon_T = _total_strain[_qp];
  epsilon_e = _elastic_strain[_qp];
  epsilon_p = _plastic_strain[_qp]; // from e_T = e_e + e_p
  epsilon_m = _mechanical_strain[_qp];

  epsilon_T_dot = (_total_strain[_qp] - _total_strain_old[_qp]) / _dt;
  epsilon_e_dot = (_elastic_strain[_qp] - _elastic_strain_old[_qp]) / _dt;
  epsilon_p_dot = (_plastic_strain[_qp] - _plastic_strain_old[_qp]) / _dt; // by conventional derivative additive property
  epsilon_m_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  
  Real epsilon_T_trace_dot, epsilon_e_trace_dot, epsilon_p_trace_dot, epsilon_m_trace_dot;
  epsilon_T_trace_dot = (_total_strain[_qp].trace() - _total_strain_old[_qp].trace()) / _dt;
  epsilon_e_trace_dot =	(_elastic_strain[_qp].trace() - _elastic_strain_old[_qp].trace()) / _dt;
  epsilon_p_trace_dot = (_plastic_strain[_qp].trace() - _plastic_strain_old[_qp].trace()) / _dt;
  epsilon_m_trace_dot = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;

  // compute artificial viscosity stress and heat source terms
  Real q_tot;
  q_tot = 0.0;

  Real stress_av, q_av, J;
  // stress_av = _density[_qp] * _C0 * std::pow(_h_max, 2) * (std::pow(epsilon_T_dot.trace(), 2) - _C1 * _ss_prop[_qp] * epsilon_T_dot.trace());
  J = 1 + _mechanical_strain[_qp].trace();
  stress_av = ( _C0 * _density[_qp] * epsilon_m_trace_dot * std::abs(epsilon_m_trace_dot) * std::pow(_h_max,2.0) / std::pow(J,2.0) ) + ( _C1 * _density[_qp] * _ss_prop[_qp] * epsilon_m_trace_dot * _h_max / J );
  q_av = _beta_av * stress_av * epsilon_m_dot.trace();
  
  if (J < 1) {
      q_tot += std::abs(q_av);  // to ensure that only when in compression the av term contributes to heat generation
  }
  
  // compute plastic flow stress and heat source
  //Real q_p;
  RankTwoTensor stress_e, q_p;
  stress_e = _elasticity_tensor[_qp] * (_mechanical_strain[_qp] - _plastic_strain[_qp]); // elastic-only stress
  q_p = _beta_p * stress_e * epsilon_p_dot;
  q_tot += std::abs(q_p.doubleContraction(I2));

  // compute thermomechanical coupled term
  Real temp, e0, alpha;
  temp = _u[_qp]; // temperature variable

  q_eos = - temp * _gamma * _density[_qp] * _specific_heat[_qp] * epsilon_m_dot.trace();
  q_tot += q_eos;
  
  return - std::abs(q_eos) * _phi[_j][_qp] * _test[_i][_qp];
}
