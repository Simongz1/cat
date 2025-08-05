#include "ADPressureHS.h"

registerMooseObject("catApp", ADPressureHS);

InputParameters
ADPressureHS::validParams()
{
  InputParameters params = ADMatHeatSource::validParams();
  params.addClassDescription("compute elastic thermomechanical coupling heat source term");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredCoupledVar("Yfinal", "reacted mass fraction");
  params.addRequiredParam<Real>("beta_av", "beta_av");
  return params;
}

ADPressureHS::ADPressureHS(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    _T(coupledValue("temperature")),
    _T_var_number(coupled("temperature")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _Yfinal(coupledValue("Yfinal")),
    _pressure_mg(getADMaterialProperty<Real>("pressure_mg")),
    _pressure_JWL(getADMaterialProperty<Real>("pressure_JWL")),
    _beta_av(getParam<Real>("beta_av")),
    _pressure_av(getADMaterialProperty<Real>("pressure_av")),
    _dP_dT(getADMaterialProperty<Real>("dP_dT")),
    _Ee_dot(getMaterialProperty<RankTwoTensor>("Ee_dot")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient"))
{}

ADReal
ADPressureHS::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  ///full pressure

  ADReal res_pressure;
  res_pressure = - _u[_qp] * _dP_dT[_qp] * _Ee_dot[_qp].trace();

  ADReal res_av;
  res_av = _beta_av * _pressure_av[_qp] * _Ee_dot[_qp].trace();

  ADReal res_tot = res_pressure + res_av;
  return - res_tot * _test[_i][_qp];
}