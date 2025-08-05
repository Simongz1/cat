#include "EosHeatSource.h"

registerMooseObject("beaverApp",EosHeatSource);

//template <>
InputParameters
//validParams<ThermalExpansionHeatSourceFiniteStrainMieGruneisenNew>()
EosHeatSource::validParams()
{
  //InputParameters params = validParams<HeatSource>();
  InputParameters params = HeatSource::validParams();
  params.addClassDescription("Thermal expansion heat source kernel"
                             "generic kernel for finite strain"
                             "for Mie Gruneisen equation of state"
                             "(Menon, 2014) (Zhang, 2011)");
  // params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  // params.addRequiredParam<Real>("s_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  // params.addCoupledVar("c","Phase field damage variable: Used to indicate calculation of Off Diagonal Jacobian term");
  params.addRequiredParam<Real>("h_max","Maximum element size");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("beta_av", "artificial viscosity beta parameter");
  params.addRequiredParam<Real>("beta_p", "plastic flow beta parameter");
  params.addRequiredParam<Real>("ss", "sound speed");
  return params;
}

EosHeatSource::EosHeatSource(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _G_Gruneisen(getParam<Real>("G_Gruneisen")), // Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density_name")),
    _deformation_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient
    _deformation_gradient_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient, previous timestep
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")), //elasticity tensor
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _total_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain")),
    _total_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "total_strain")),
    _elastic_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "total_strain")),
    // _PK1(declareProperty<RankTwoTensor>("PK1")),
    // _PK2(declareProperty<RankTwoTensor>("PK2")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _beta_av(getParam<Real>("beta_av")),
    _beta_p(getParam<Real>("beta_p")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _h_max(getParam<Real>("h_max")),
    _ss(getParam<Real>("ss"))

    // declare material properties to plot PK stresses and heat sources separately

    // stresses
    //_PK1(declareProperty<RankTwoTensor>("PK1")),
    //_PK2(declareProperty<RankTwoTensor>("PK2")),

    // heat sources

    //_q_eos(declareProperty<Real>("q_eos")),
    //_q_cpl(declareProperty<Real>("q_cpl")),
    //_q_av(declareProperty<Real>("q_av")),
    //_q_p(declareProperty<Real>("q_p")),
    //_q_tot(declareProperty<Real>("q_tot"))

{

}

Real
EosHeatSource::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  Real heat_source, J, J_dot;
  Real Kb = _Bulk_Modulus_Ref;
  Real thermal_expansion_coeff, q_eos; // thermal expansion coefficient depends on Gruneisen parameter, bulk modulus and sound speed
  RankTwoTensor Ce, Ce_old, Ce_inv, strain_rate;
  RankTwoTensor F, F_old, F_dot, F_inv;

  // calculate deformation gradient
  F = _deformation_gradient[_qp];
  F_old = _deformation_gradient_old[_qp];
  F_dot = (F - F_old) / _dt;

  J = F.det();
  J_dot = F_dot.det();

  F_inv = F.inverse();

  // preliminary calculation of deformation-related tensors
  Ce = F.transpose() * F; // Ce = Cauchy tensor
  Ce_old = F_old.transpose() * F_old; // Ce old = Cauchy tensor at previous time step
  Ce_inv = Ce.inverse();

  strain_rate = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;

  // decomposition of EOS heat source into eso and cpl components

  Real temp;
  temp = _u[_qp];
  q_eos = - temp * _G_Gruneisen * _density[_qp] * _specific_heat[_qp] * (Ce_inv * strain_rate).trace();
  //_q_eos[_qp] = q_eos;

  thermal_expansion_coeff = _G_Gruneisen * _density[_qp] * _specific_heat[_qp] / Kb;

  // compute coupling component

  RankTwoTensor thermal_coupling_tensor;
  Real q_cpl, q_tot;

  thermal_coupling_tensor = _elasticity_tensor[_qp] * I2;
  thermal_coupling_tensor = strain_rate * thermal_coupling_tensor;
  q_cpl = thermal_expansion_coeff * temp
              * std::exp((2.0/3.0) * thermal_expansion_coeff * (temp - _reference_temperature))
              * (Kb * std::pow(J , 2.0/3.0) * (Ce_inv * strain_rate).trace() - (1.0/3.0) * thermal_coupling_tensor.trace());

  //_q_cpl[_qp] = q_cpl;       
  // compute elastic and plastic lagrangian strains
  
  RankTwoTensor Ee, Ee_dot, Ep, Ep_dot, Et, Et_dot;

  Ee = _elastic_strain[_qp];
  Ee_dot = (_elastic_strain[_qp] - _elastic_strain_old[_qp]) / _dt;

  Et = _total_strain[_qp];
  Et_dot = (_total_strain[_qp] - _total_strain_old[_qp]) / _dt;

  Ep = Et - Ee;
  Ep_dot = Et_dot - Ee_dot;

  // compute PK1 and PK2 stresses

  RankTwoTensor PK1, PK2, Cauchy_stress;

  Cauchy_stress = _stress[_qp];
  PK1 = J * Cauchy_stress * F_inv;
  PK2 = J * F_inv * Cauchy_stress * (F_inv.transpose()); // this stress already contains eos + av + pf as it comes from eos stress material

  // append into the material properties created

  //_PK1[_qp] = PK1;
  //_PK2[_qp] = PK2;
  
  // compute av contribution

  Real q_av;
  RankTwoTensor stress_av;

  stress_av = (_C0 * _density[_qp] * J_dot * std::abs(J_dot) * (1.0 / std::pow(J, 3)) * std::pow(_h_max, 2)) * I2;
  stress_av += (_C1 * _density[_qp] * _ss * (J_dot / std::pow(J, 2)) * _h_max) * I2;

  q_av = std::abs(_beta_av * stress_av.doubleContraction(Ee_dot));
  //_q_av[_qp] = q_av;

  // compute plastic flow heat

  Real q_p;

  q_p = std::abs(_beta_p * PK2.doubleContraction(Ep_dot));
  //_q_p[_qp] = q_p;
  
  q_tot = q_eos + q_cpl + q_av + q_p;
  //_q_tot[_qp] = q_tot;

  return - q_tot * _test[_i][_qp];
}

Real
EosHeatSource::computeQpJacobian()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  Real heat_source, J, J_dot;
  Real Kb = _Bulk_Modulus_Ref;
  Real thermal_expansion_coeff, q_eos; // thermal expansion coefficient depends on Gruneisen parameter, bulk modulus and sound speed
  RankTwoTensor Ce, Ce_old, Ce_inv, strain_rate;
  RankTwoTensor F, F_old, F_dot, F_inv;

  // calculate elastic deformation gradient
  F = _deformation_gradient[_qp];
  F_old = _deformation_gradient_old[_qp];
  F_dot = (F - F_old) / _dt;

  J = F.det();
  J_dot = F_dot.det();

  F_inv = F.inverse();

  // preliminary calculation of deformation-related tensors
  Ce = F.transpose() * F; // Ce = Cauchy tensor
  Ce_old = F_old.transpose() * F_old; // Ce old = Cauchy tensor at previous time step
  Ce_inv = Ce.inverse();

  strain_rate = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;

  // decomposition of EOS heat source into eso and cpl components

  Real temp;
  temp = _u[_qp];
  q_eos = - temp * _G_Gruneisen * _density[_qp] * _specific_heat[_qp] * (Ce_inv * strain_rate).trace();
  //_q_eos[_qp] = q_eos;

  thermal_expansion_coeff = _G_Gruneisen * _density[_qp] * _specific_heat[_qp] / Kb;

  // compute coupling component

  RankTwoTensor thermal_coupling_tensor;
  Real q_cpl, q_tot;

  thermal_coupling_tensor = _elasticity_tensor[_qp] * I2;
  thermal_coupling_tensor = strain_rate * thermal_coupling_tensor;
  q_cpl = - thermal_expansion_coeff * temp
              * std::exp((2.0/3.0) * thermal_expansion_coeff * (temp - _reference_temperature))
              * (Kb * std::pow(J , 2.0/3.0) * (Ce_inv * strain_rate).trace() - (1.0/3.0) * thermal_coupling_tensor.trace());

  //_q_cpl[_qp] = q_cpl;

  // compute elastic and plastic strains
  
  RankTwoTensor Ee, Ee_dot, Ep, Ep_dot, Et, Et_dot;

  Ee = _elastic_strain[_qp];
  Ee_dot = (_elastic_strain[_qp] - _elastic_strain_old[_qp]) / _dt;

  Et = _total_strain[_qp];
  Et_dot = (_total_strain[_qp] - _total_strain_old[_qp]) / _dt;

  Ep = Et - Ee;
  Ep_dot = Et_dot - Ee_dot;

  // compute PK1 and PK2 stresses

  RankTwoTensor PK1, PK2, Cauchy_stress;

  Cauchy_stress = _stress[_qp];
  PK1 = J * Cauchy_stress * F_inv;
  PK2 = J * F_inv * Cauchy_stress * (F_inv.transpose()); // this stress already contains eos + av + pf as it comes from eos stress material

  // append into the material properties created

  //_PK1[_qp] = PK1;
  //_PK2[_qp] = PK2;

  // compute av contribution

  Real q_av;
  RankTwoTensor stress_av, PK2_c;

  stress_av = (_C0 * _density[_qp] * J_dot * std::abs(J_dot) * (1.0 / std::pow(J, 3)) * std::pow(_h_max, 2)) * I2;
  stress_av += (_C1 * _density[_qp] * _ss * (J_dot / std::pow(J, 2)) * _h_max) * I2;

  q_av = _beta_av * stress_av.doubleContraction(Ee_dot);
  //_q_av[_qp] = q_av;

  PK2_c = PK2 - stress_av * I2; // decomposition of PK2 stress into conservative and dissipative part (av)

  // compute plastic flow heat

  Real q_p;

  q_p = _beta_p * PK2_c.doubleContraction(Ep_dot);
  //_q_p[_qp] = q_p;
  
  q_tot = q_eos + q_cpl + q_av + q_p;
  //_q_tot[_qp] = q_tot;

  return - std::abs(_G_Gruneisen * _density[_qp] * _specific_heat[_qp] * (Ce_inv * strain_rate).trace()) * _phi[_j][_qp] * _test[_i][_qp];
}
