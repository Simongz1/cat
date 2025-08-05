#include "MGHeatSourceV2.h"

registerMooseObject("beaverApp", MGHeatSourceV2);

InputParameters
MGHeatSourceV2::validParams()
{
  InputParameters params = HeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("Heat Source contribution of"
                              "Mie Gruneisen Thermo Mechanical"
                              "Coupling.");
  params.addRequiredParam<Real>("gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  params.addRequiredParam<Real>("beta_av", "artificial viscosity beta parameter");
  return params;
}

MGHeatSourceV2::MGHeatSourceV2(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _gamma(getParam<Real>("gamma")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density_name")),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")), //elasticity tensor
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>("mechanical_strain")),
    _elastic_strain(getMaterialProperty<RankTwoTensor>("elastic_strain")),
    _elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("elastic_strain")),
    _beta_av(getParam<Real>("beta_av")),
    _stress_av(getMaterialProperty<Real>("stress_av"))
{}

Real
MGHeatSourceV2::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  Real heat_source;
  Real thermal_expansion_coeff; // thermal expansion coefficient depends on Gruneisen parameter, bulk modulus and sound speed
  RankTwoTensor epsilon_e, epsilon_m, epsilon_m_dot;

  // compute total, elastic and plastic infinitesimal strains with time derivatives
  epsilon_e = _elastic_strain[_qp];
  epsilon_m = _mechanical_strain[_qp];

  RankTwoTensor epsilon_e_dot = (_elastic_strain[_qp] - _elastic_strain_old[_qp]) / _dt;
  epsilon_m_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;

  // initialize the heat source value
  Real q_tot;
  q_tot = 0.0;

  // compute the artificial viscosity contribution
  Real q_av;
  Real stress_av = _stress_av[_qp]; //extract the QP value of the AV stress
  q_av = _beta_av * stress_av * epsilon_m_dot.trace();
  q_tot += std::abs(q_av); //append the value of the heat source due to AV to the total HS value

  // compute thermomechanical coupling term
  Real e0, alpha;
  Real q_eos = - _u[_qp] * _gamma * _density[_qp] * _specific_heat[_qp] * epsilon_e_dot.trace(); //compute EOS thermomechanical coupling
  q_tot += std::abs(q_eos);

  return - q_tot * _test[_i][_qp]; //residual QP value
}
