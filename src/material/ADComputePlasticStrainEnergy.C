#include "ADComputePlasticStrainEnergy.h"

registerMooseObject("mlApp", ADComputePlasticStrainEnergy);

InputParameters
ADComputePlasticStrainEnergy::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADMaterial>::validParams();
  //params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("plastic strain rate constitutive model");
  params.addRequiredParam<Real>("N_p", "N_p");
  params.addRequiredParam<Real>("mu_p", "mu_p");
  //variables

  //parameters
  
  return params;
}

ADComputePlasticStrainEnergy::ADComputePlasticStrainEnergy(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ADMaterial>(parameters),
    //ADSingleVariableReturnMappingSolution(parameters),
    _epsilon_p(getADMaterialProperty<RankTwoTensor>("epsilon_p")),
    _N_p(getParam<Real>("N_p")),
    _mu_p(getParam<Real>("mu_p")),
    _exp2ep(declareADProperty<RankTwoTensor>("exp2ep")),
    _Wp(declareADProperty<Real>("Wp")),
    _sp(declareADProperty<RankTwoTensor>("sp"))
{}

void
ADComputePlasticStrainEnergy::initialSetup()
{}

void
ADComputePlasticStrainEnergy::initQpStatefulProperties()
{
  ADMaterial::initQpStatefulProperties();
}

void
ADComputePlasticStrainEnergy::computeQpProperties()
{
  //compute spectral decomposition of strain
  ADRankTwoTensor Q;
  std::vector<ADReal> lam(3);
  
  //get decomposition of 2 * epsilon_p
  ADRankTwoTensor eps2 = 2. * _epsilon_p[_qp];
  eps2.symmetricEigenvaluesEigenvectors(lam, Q);

  //assemble exponential matrix
  ADRankTwoTensor expD;
  expD.zero();

  for(unsigned int i = 0; i < 3; ++i){
    expD(i,i) = MetaPhysicL::exp(lam[i]);
  }

  //assemble the total matrix for exp(2*eps)
  ADRankTwoTensor exp2ep = Q * expD * Q.transpose();

  //store
  _exp2ep[_qp] = exp2ep;

  //now that we have this, we can proceed to compute the energy and the derivative

  //ENERGY COMPUTATION
  ADReal lambda_p = MetaPhysicL::sqrt(1. / 3. * _exp2ep[_qp].trace());
  ADReal lambda_rp = lambda_p / std::sqrt(_N_p);

  ADReal Wp_term1 = lambda_rp * computeInverseLang(lambda_rp);
  ADReal Wp_term2 = MetaPhysicL::log(computeInverseLang(lambda_rp) / MetaPhysicL::sinh(computeInverseLang(lambda_rp)));
  _Wp[_qp] = _mu_p * _N_p * (Wp_term1 + Wp_term2);

  //compute the derivative of the plastic energy with respect to plastic strain
  _sp[_qp] = computedWdInvLang(lambda_rp) * computeInverseLangDerivative(lambda_rp) * computedLambdadep();
}

ADReal
ADComputePlasticStrainEnergy::computeInverseLang(const ADReal &x)
{
  return x * (3. - MetaPhysicL::pow(x, 2.)) / (1. - MetaPhysicL::pow(x, 2.));
}

ADReal
ADComputePlasticStrainEnergy::computeInverseLangDerivative(const ADReal &x)
{
  ADReal num = (3. - 3. * MetaPhysicL::pow(x, 2.)) * (1. - MetaPhysicL::pow(x, 2.));
  num += x * (3. - MetaPhysicL::pow(x, 2.));

  ADReal den = 1. - MetaPhysicL::pow(x, 2.);
  return num / den;
}

ADReal
ADComputePlasticStrainEnergy::computedWdInvLang(const ADReal &x)
{
  ADReal term1 = _mu_p * _N_p * x;
  ADReal term2 = MetaPhysicL::sinh(computeInverseLang(x)) / computeInverseLang(x);
  ADReal term3 = - 1. * MetaPhysicL::cosh(computeInverseLang(x));
  return term1 + term2 + term3;
}

ADRankTwoTensor
ADComputePlasticStrainEnergy::computedLambdadep()
{
  ADRankTwoTensor I2;
  I2.setToIdentity();

  ADRankTwoTensor res = 1. / (6. * std::sqrt(_N_p)) * 2. * _exp2ep[_qp];
  res *= 1. / MetaPhysicL::sqrt(1. / 3. * _exp2ep[_qp].trace());
  return res;
}