#include "ADVisHSLIPIT.h"

registerMooseObject("mistApp", ADVisHSLIPIT);

InputParameters
ADVisHSLIPIT::validParams()
{
  InputParameters params = ADMatHeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute heating from viscous dissipation");
  params.addRequiredParam<Real>("beta_p", "beta_p");
  params.addRequiredParam<Real>("beta_comp", "beta_compression");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature");
  return params;
}

ADVisHSLIPIT::ADVisHSLIPIT(const InputParameters & parameters)
  : ADMatHeatSource(parameters),
    //retrieve Ydots
    _cauchy_stress(getMaterialProperty<RankTwoTensor>("cauchy_stress")),
    _beta_p(getParam<Real>("beta_p")),
    _beta_comp(getParam<Real>("beta_comp")),
    _Ep_dot(getMaterialProperty<RankTwoTensor>("Ep_dot")),
    _Ee_dot(getMaterialProperty<RankTwoTensor>("Ee_dot")),
    _sigma(getMaterialProperty<RankTwoTensor>("sigma")),
    _F(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _ep_rate(getMaterialProperty<Real>("ep_rate")),
    _alpha(getADMaterialProperty<Real>("alpha")),
    _Fe(getMaterialProperty<RankTwoTensor>("Fe")),
    _Fe_old(getMaterialPropertyOld<RankTwoTensor>("Fe")),
    _Tref(getParam<Real>("reference_temperature")),
    _Cijkl(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _S(getMaterialProperty<RankTwoTensor>("S")),
    _HS_plastic(getMaterialProperty<Real>("HS_plastic")),
    _C_computed(getMaterialProperty<RankTwoTensor>("C_computed"))
{}

ADReal
ADVisHSLIPIT::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  ADReal q_plastic = _beta_p * _HS_plastic[_qp];

  //include heating due to compression
  Real K = ElasticityTensorTools::getIsotropicBulkModulus(_Cijkl[_qp]);
  ADReal q_compression = - _beta_comp * K * _alpha[_qp] * (_u[_qp]) * (_C_computed[_qp].inverse().doubleContraction(_Ee_dot[_qp]));

  return - (q_plastic + std::max(q_compression, 0.)) * _test[_i][_qp];
}