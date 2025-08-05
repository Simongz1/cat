#include "ComputeLagrangianAVStress.h"

registerMooseObject("TensorMechanicsApp", ComputeLagrangianAVStress);

InputParameters
ComputeLagrangianAVStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes Lagrangian form of artificial viscosity");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredParam<Real>("element_size", "element_size");
  return params;
}

ComputeLagrangianAVStress::ComputeLagrangianAVStress(
    const InputParameters & parameters)
  : Material(parameters),
    _rho(getMaterialProperty<Real>("density")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("element_size")),
    _pressure_av(declareProperty<RankTwoTensor>("pressure_av")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor"))
{
}

void
ComputeLagrangianAVStress::computeQpProperties()
{
  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour
  const Real K = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  Real K0 = K;
  Real ss = std::sqrt(K / _rho[_qp]);
	
  //Compute artificial viscosity term
  Real P_av;
  Real Je;
  Real Je_dot;
  Je = _deformation_gradient[_qp].det();
  Je_dot = ((_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det()) / _dt);

  P_av = _C0 * _rho[_qp] * (Je_dot * std::abs(Je_dot) / std::pow(Je, 2.0)) * std::pow(_Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (Je_dot / Je) * _Le;
  _pressure_av[_qp] = P_av * I2;
}
