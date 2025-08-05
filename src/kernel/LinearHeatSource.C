// HEAT SOURCE DUE TO PLASTIC STRAIN J2 PLASTICITY and EOS contribution

#include "LinearHeatSource.h"
#include <tuple>

registerMooseObject("beaverApp", LinearHeatSource);

InputParameters
LinearHeatSource::validParams()
{
  InputParameters params = HeatSource::validParams();
  params.addClassDescription("Linear elastic thermal coupling");
  params.addRequiredCoupledVar("temperature", "temperature field");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  params.addRequiredParam<Real>("alpha", "thermal expansion coeff");
  return params;
}

LinearHeatSource::LinearHeatSource(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),

    _temperature(coupledValue("temperature")),
    // _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density")),

    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")), //elasticity tensor
    _alpha(getParam<Real>("alpha")),
        
    _total_strain_old(getMaterialPropertyOldByName<RankTwoTensor>("total_strain")),
    _total_strain(getMaterialPropertyByName<RankTwoTensor>("total_strain")),
    
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>("mechanical_strain")),
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress"))
{
}

Real
LinearHeatSource::computeQpResidual()
{
  Real heat_source_density, q_source;
  Real temperature, ref;
  RankTwoTensor C, C_inv, stress, total_strain, strain_rate;
  Real Res1, Res2;

  RankTwoTensor I2(RankTwoTensor::initIdentity);

  // K = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  // J = det(F)

  strain_rate = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;

  q_source = - _alpha * _temperature[_qp] * strain_rate.doubleContraction(I2);

  // compute second term (EOS + CPL terms)

  // additional term for virtual viscosity

  //delta_strain = _effective_plastic_strain[_qp] - _effective_plastic_strain_old[_qp];
  //temperature_new = temperature + ((_eta) / (_alpha)) * delta_strain;

  //delta_temp_ref = temperature - _reference_temperature;
  //delta_temp = temperature_new - temperature;

  // thermal_strain = _alpha * delta_temp;
  //heat_source_density = _density[_qp] * _specific_heat[_qp] * delta_temp * vol_strain;
  //effective_ps = _effective_plastic_strain[_qp];
  // effective_ps = _effective_ps;

  Res1 = - q_source * _test[_i][_qp];
  return Res1;
}

Real
LinearHeatSource::computeQpJacobian()
{
  Real heat_source_density, q_source;
  Real temperature, ref;
  RankTwoTensor C, C_inv, stress, total_strain, strain_rate;
  Real Res1, Res2;

  RankTwoTensor I2(RankTwoTensor::initIdentity);

  // K = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  // J = det(F)

  strain_rate = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;

  q_source = - _alpha * _temperature[_qp] * strain_rate.doubleContraction(I2);

  // compute second term (EOS + CPL terms)

  // additional term for virtual viscosity

  //delta_strain = _effective_plastic_strain[_qp] - _effective_plastic_strain_old[_qp];
  //temperature_new = temperature + ((_eta) / (_alpha)) * delta_strain;

  //delta_temp_ref = temperature - _reference_temperature;
  //delta_temp = temperature_new - temperature;

  // thermal_strain = _alpha * delta_temp;
  //heat_source_density = _density[_qp] * _specific_heat[_qp] * delta_temp * vol_strain;
  //effective_ps = _effective_plastic_strain[_qp];
  // effective_ps = _effective_ps;

  Res2 = - q_source * _phi[_j][_qp] * _test[_i][_qp];;
  return Res2;
}