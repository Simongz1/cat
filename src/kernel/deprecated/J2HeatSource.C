// HEAT SOURCE DUE TO PLASTIC STRAIN J2 PLASTICITY and EOS contribution

#include "J2HeatSource.h"
#include <tuple>

registerMooseObject("beaverApp", J2HeatSource);

InputParameters
J2HeatSource::validParams()
{
  InputParameters params = HeatSource::validParams();
  params.addClassDescription("J2 plasticity heat source kernel");
  params.addRequiredCoupledVar("temperature", "temperature field");
  params.addRequiredCoupledVar("effective_ps", "effective plastic strain");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  params.addRequiredParam<Real>("alpha", "thermal expansion coeff");
  // params.addRequiredParam<MaterialPropertyName>("effective_plastic_strain", "the effective plastic strain from J2 plasticity");
  params.addRequiredParam<Real>("eta", "thermal conversion parameter");
  params.addRequiredParam<Real>("beta_p", "beta coefficient for plastic contribution");
  params.addRequiredParam<Real>("beta_vis", "beta coefficient for viscosity contribution");
  params.addRequiredParam<Real>("gamma", "gruneisen parameter");
  params.addRequiredCoupledVar("vol_strain", "vol_strain");
  params.addRequiredParam<Real>("C0", "C0");
  params.addRequiredParam<Real>("C1", "C1");
  params.addRequiredParam<Real>("h", "minimum element size");
  params.addRequiredParam<Real>("ss", "linear sound speed");
  return params;
}

J2HeatSource::J2HeatSource(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),

    _temperature(coupledValue("temperature")),
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density")),

    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")), //elasticity tensor
    _alpha(getParam<Real>("alpha")),
    _eta(getParam<Real>("eta")),
    _effective_plastic_strain_old(getMaterialPropertyOldByName<Real>("effective_plastic_strain")),
    _effective_plastic_strain(getMaterialPropertyByName<Real>("effective_plastic_strain")),
    _effective_ps(coupledValue("effective_ps")),
    _vol_strain(coupledValue("vol_strain")),
    _pk1_stress(getMaterialPropertyByName<RankTwoTensor>("pk1_stress")),
    // _pk2_stress(getMaterialPropertyByName<RankTwoTensor>("pk2_stress")),
    _beta_p(getParam<Real>("beta_p")),
    _total_strain_old(getMaterialPropertyOldByName<RankTwoTensor>("total_strain")),
    _total_strain(getMaterialPropertyByName<RankTwoTensor>("total_strain")),
    _gamma(getParam<Real>("gamma")),
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>("mechanical_strain")),
    _deformation_gradient(getMaterialPropertyByName<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOldByName<RankTwoTensor>("deformation_gradient")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _h(getParam<Real>("h")),
    _ss(getParam<Real>("ss")),
    _beta_vis(getParam<Real>("beta_vis"))
    // _effective_ps(declareProperty<Real>("effective_ps"))
{

}

Real
J2HeatSource::computeQpResidual()
{
  Real heat_source_density, K;
  Real temperature, delta_temp, delta_temp_ref, temperature_new, ref;
  Real delta_strain, thermal_strain;
  Real effective_ps, vol_strain;
  Real Res1, Res2;
  RankTwoTensor pk1, pk2, E_p_dot, E_e_dot, F, F_dot, F_inv, S_vis, E_dot, ce, ce_inv, cc, cc_inv;
  RankFourTensor C, C_inv;
  Real beta_p, J, q_pk2, q_cpl, q_eos, q_tot; //beta coefficient for the contribution of plastic stress to heat

  RankTwoTensor I2(RankTwoTensor::initIdentity);

  F = _deformation_gradient[_qp];
  F_inv = F.inverse();
  pk1 = _pk1_stress[_qp];
  pk2 = F_inv * pk1;
  K = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  J = F.det();
  C = _elasticity_tensor[_qp];
  ce = (1 / J) * pk1 * F; // elastic cauchy strain tensor
  cc = (F.transpose()) * (F); // Cauchy strain tensor
  cc_inv = cc.inverse(); // cc inverse tensor
  ce_inv = ce.inverse();
  C_inv = C.inverse();
  vol_strain = _vol_strain[_qp];
  beta_p = _beta_p;
  E_p_dot = (_total_strain[_qp] - _total_strain_old[_qp]) / _dt; // time derivative of plastic lagrangian strain
  E_e_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt; // time derivative of elastic lagrangian strain
  E_dot = E_p_dot + E_e_dot; //total lagrangian strain time derivative
  temperature = _temperature[_qp];
  ref = _reference_temperature;

  // compute first term beta * PK2 * E_p_dot

  // pk2 = _pk2_stress[_qp];
  q_pk2 = beta_p * pk2.doubleContraction(E_p_dot);

  // compute second term (EOS + CPL terms)

  q_eos = - _gamma * _density[_qp] * _specific_heat[_qp] * temperature * (cc_inv * E_e_dot).trace();
  q_cpl = - (_alpha * temperature / 3) * std::exp((2.0 / 3.0) * _alpha * (temperature - ref)) * (E_e_dot.doubleContraction(C * I2)) + 
	K * _alpha * temperature * std::pow(J, 2.0 / 3.0) * std::exp(((2.0 / 3.0) * _alpha * (temperature - ref))) * (cc_inv * E_e_dot).trace();

  q_tot = q_pk2 + q_eos + q_cpl;
  
  // additional virtual viscosity term
  
  F_dot = (_deformation_gradient[_qp] - _deformation_gradient_old[_qp]) / _dt;
  Real J_dot;
  J_dot = F_dot.det();
  S_vis = _C0 * _density[_qp] * ((J_dot * std::abs(J_dot)) / J * J) * (_h * _h) * I2;
  S_vis += _C1 * _density[_qp] * _ss * (J_dot / J) * _h * I2;

  // compute viscosity heat contribution
  
  Real q_vis;
  q_vis = _beta_vis * S_vis.doubleContraction(E_dot * I2);
  
  // update heat source
  
  q_tot += q_vis;

  Res1 = q_tot * _test[_i][_qp];
  Res2 = effective_ps * _test[_i][_qp];

  return Res1;
}
