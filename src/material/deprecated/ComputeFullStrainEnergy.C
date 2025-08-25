#include "ComputeFullStrainEnergy.h"

registerMooseObject("catApp", ComputeFullStrainEnergy);

InputParameters
ComputeFullStrainEnergy::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes the strain energy density using full elasticity tensor");
  params.addRequiredParam<Real>("T0", "stress free temperature");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("poisson_ratio", "poisson_ratio");
  params.addRequiredParam<Real>("beta", "beta parameter for blatz - ko");
  return params;
}

ComputeFullStrainEnergy::ComputeFullStrainEnergy(const InputParameters & parameters)
  : Material(parameters),
    _alpha(getMaterialProperty<Real>("alpha")),
    _W_tot(declareProperty<Real>("W_tot")), //total strain energy density
    _W_el(declareProperty<Real>("W_el")), //elastic strain energy density
    _W_th(declareProperty<Real>("W_th")), //thermal strain energy density
    _W_bk(declareProperty<Real>("W_bk")), //blatz-ko strain energy density

    _dW_bk_dIv1(declareProperty<Real>("dW_bk_dIv1")),
    _dW_bk_dIv2(declareProperty<Real>("dW_bk_dIv2")),
    _dW_bk_dIv3(declareProperty<Real>("dW_bk_dIv3")),

    //output invariants
    _Iv1(declareProperty<Real>("Iv1")),
    _Iv2(declareProperty<Real>("Iv2")),
    _Iv3(declareProperty<Real>("Iv3")),

    //retrieve required tensors
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _T0(getParam<Real>("T0")),
    _temperature(coupledValue("temperature")),
    _nu(getParam<Real>("poisson_ratio")),
    _beta(getParam<Real>("beta")),
    _cauchy_bk(declareProperty<RankTwoTensor>("cauchy_bk")),
    _lagrangian_strain(declareProperty<RankTwoTensor>("lagrangian_strain")),
    _lagrangian_strain_rate(declareProperty<RankTwoTensor>("lagrangian_strain_rate"))
{}

void
ComputeFullStrainEnergy::computeQpProperties()
{
  Real W_tot, W_el, W_th, W_frac;
  RankTwoTensor C, Cinv, F, E, B;

  //compute lagrangian strain
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  F = _deformation_gradient[_qp];
  RankTwoTensor F_dot = (_deformation_gradient[_qp] - _deformation_gradient_old[_qp]) / _dt;
  C = F.transpose() * F;
  B = F * F.transpose();
  E = (1.0 / 2.0) * (C - I2);
  RankTwoTensor E_dot = 0.5 * ((F_dot.transpose() * F) + (F.transpose() * F_dot));

  _lagrangian_strain[_qp] = E;
  _lagrangian_strain_rate[_qp] = E_dot;

  const Real mu = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  //compute energies
  RankTwoTensor W_el_right;
  W_el_right = (_elasticity_tensor[_qp] * I2).doubleContraction(E); //right part of contraction
  W_el = (1.0 / 2.0) * (E.doubleContraction(W_el_right));

  RankTwoTensor W_th_right, alpha_tensor;
  alpha_tensor = _alpha[_qp] * I2 * (_temperature[_qp] - _T0);
  W_th_right = (_elasticity_tensor[_qp] * I2).doubleContraction(alpha_tensor);
  W_th = E.doubleContraction(W_th_right);

  W_tot = W_el - W_th;

  //write into properties

  _W_tot[_qp] = W_tot;
  _W_el[_qp] = W_el;
  _W_th[_qp] = W_th;

  //compute blatz ko strain energy and derivatives

  Real W_bk;
  Real Iv1, Iv2, Iv3;

  //define invariants and write into material properties to be called by stress material
  Iv1 = C.trace();
  _Iv1[_qp] = Iv1;

  Iv2 = 0.5 * (std::pow(C.trace(), 2) - (C * C).trace());
  _Iv2[_qp] = Iv2;

  Iv3 = C.det();
  _Iv3[_qp] = Iv3;

  Real a = - std::pow((1.0 - 2.0 * _nu), - 1.0);

  //W_bk = (mu / 2.0) * (Iv1 + (1.0 / a) * std::pow(Iv3, - a) - 3.0);
  //W_bk += (mu / 2.0) * (1.0 - _beta) * ((Iv2 / Iv3) + (1.0 / a) * (std::pow(Iv3, a) - 1.0) - 3.0);
  W_bk = mu * ((Iv2 / Iv3) + 2.0 * std::sqrt(Iv3) - 5);
  _W_bk[_qp] = W_bk;

  //compute derivatives
  Real dW_bk_dIv1, dW_bk_dIv2, dW_bk_dIv3;

  //dW_bk_dIv1 = mu / 2.0;
  //dW_bk_dIv2 = (mu / 2.0) * (1.0 - _beta) / Iv3;
  //dW_bk_dIv3 = (mu / 2.0) * (- std::pow(Iv3, - (a + 1.0)));
  //dW_bk_dIv3 += (mu / 2.0) * (1.0 - _beta) * (- (Iv2 / std::pow(Iv3, 2.0)) + std::pow(Iv3, a - 1.0));
  dW_bk_dIv1 = 0.0;
  dW_bk_dIv2 = 1.0 / Iv3;
  dW_bk_dIv3 = - (Iv2 / Iv3) + (1.0 / std::sqrt(Iv3)); 
  
  _dW_bk_dIv1[_qp] = dW_bk_dIv1;
  _dW_bk_dIv2[_qp] = dW_bk_dIv2;
  _dW_bk_dIv3[_qp] = dW_bk_dIv3;

  //provitional computation of cauchy stress from blatz - ko strain energy
  _cauchy_bk[_qp] = (2.0 / std::sqrt(Iv3)) * ((Iv3 * dW_bk_dIv3 * I2) + ((dW_bk_dIv1 + (Iv1 * dW_bk_dIv2)) * B) - dW_bk_dIv2 * (B * B));
}
