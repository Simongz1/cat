#include "LIPITHS.h"

registerMooseObject("mlApp", LIPITHS);

InputParameters
LIPITHS::validParams()
{
  InputParameters params = ADKernel::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute heating from viscous dissipation using AD");
  params.addRequiredParam<Real>("beta_p", "beta_p");
  params.addRequiredParam<Real>("beta_av", "beta_av");
  return params;
}

LIPITHS::LIPITHS(const InputParameters & parameters)
  : ADKernel(parameters),
    //retrieve Ydots
    _beta_p(getParam<Real>("beta_p")),
    _beta_av(getParam<Real>("beta_av")),
    _Ep_dot(getADMaterialProperty<RankTwoTensor>("Ep_dot")),
    _Ee_dot(getADMaterialProperty<RankTwoTensor>("Ee_dot")),
    _alpha(getADMaterialProperty<Real>("alpha")),
    _F(getADMaterialProperty<RankTwoTensor>("deformation_gradient")),

    _Cijkl(getADMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _S(getADMaterialProperty<RankTwoTensor>("S")),
    _HS_plastic(getADMaterialProperty<Real>("HS_plastic")),
    _HS_elastic(getADMaterialProperty<Real>("HS_elastic"))
{}

ADReal
LIPITHS::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  ADReal q_plastic = _beta_p * _HS_plastic[_qp];
  ADReal q_elastic = _beta_av * _HS_elastic[_qp];

  //include heating due to compression
  ADReal K = ElasticityTensorTools::getIsotropicBulkModulus(_Cijkl[_qp]);
  ADRankTwoTensor C = _F[_qp].transpose() * _F[_qp];
  ADReal q_compression = - K * _alpha[_qp] * (_u[_qp]) * (C.inverse().doubleContraction(_Ee_dot[_qp]));

  return - (q_plastic + q_elastic + std::max(q_compression, ADReal(0.))) * _test[_i][_qp];
}
