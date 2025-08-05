//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ArtVisJ2StressLIPIT.h"
#include "MathUtils.h"

registerMooseObject("catApp", ArtVisJ2StressLIPIT);

InputParameters
ArtVisJ2StressLIPIT::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ComputeLagrangianStressPK1>::validParams();
  params += SingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("The Simo-Hughes style J2 plasticity.");
  params.addParam<MaterialPropertyName>(
      "elasticity_tensor", "elasticity_tensor", "The name of the elasticity tensor.");
  params.addRequiredParam<MaterialName>("flow_stress_material",
                                        "The material defining the flow stress");
  /////////////
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredParam<Real>("element_size", "element_size");

  //fracture stuff
  params.addRequiredCoupledVar("c", "fracture variable");
  params.addRequiredParam<Real>("l", "crack length");
  params.addRequiredParam<Real>("gc", "surface energy");
  params.addRequiredCoupledVar("gcprop", "gcprop");
  params.addRequiredParam<Real>("visco", "crack viscosity parameter");
  params.addRequiredParam<Real>("kdamage", "residual stiffness");

  params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "kappa property name");
  params.addParam<MaterialPropertyName>("mobility_name", "L", "mobility property name");
  params.addParam<MaterialPropertyName>("elastic_energy_name", "elastic_energy", "elastic energy property name");
  return params;
}

ArtVisJ2StressLIPIT::ArtVisJ2StressLIPIT(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ComputeLagrangianStressPK1>(parameters),
    GuaranteeConsumer(this),
    SingleVariableReturnMappingSolution(parameters),
    _elasticity_tensor_name(_base_name + getParam<MaterialPropertyName>("elasticity_tensor")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>(_elasticity_tensor_name)),
    _F_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "deformation_gradient")),
    _ep_name(_base_name + "effective_plastic_strain"),
    _ep(declareProperty<Real>(_ep_name)),
    _ep_old(getMaterialPropertyOldByName<Real>(_ep_name)),
    _be(declareProperty<RankTwoTensor>(_base_name +
                                       "volume_preserving_elastic_left_cauchy_green_strain")),
    _be_old(getMaterialPropertyOldByName<RankTwoTensor>(
        _base_name + "volume_preserving_elastic_left_cauchy_green_strain")),
    _Np(declareProperty<RankTwoTensor>(_base_name + "flow_direction")),
    _flow_stress_material(nullptr),
    _flow_stress_name(_base_name + "flow_stress"),
    _H(getMaterialPropertyByName<Real>(_flow_stress_name)),
    _dH(getMaterialProperty<Real>("dH")),
    _d2H(getMaterialProperty<Real>("d2H")),
    /////
    _rho(getMaterialProperty<Real>("density")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("element_size")),
    _pressure_av(declareProperty<RankTwoTensor>("pressure_av")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _plastic_strain(declareProperty<RankTwoTensor>("plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("plastic_strain")),
    _plastic_strain_rate(declareProperty<RankTwoTensor>("plastic_strain_rate")),
    _pnorm(declareProperty<Real>("pnorm")),
    _hsp(declareProperty<Real>("hsp")),
    _cauchy_stress(getMaterialProperty<RankTwoTensor>("cauchy_stress")),

    //request fracture stuff
    _c(coupledValue("c")),
    _l(getParam<Real>("l")),
    _kappa(declareProperty<Real>(getParam<MaterialPropertyName>("kappa_name"))),
    _L(declareProperty<Real>(getParam<MaterialPropertyName>("mobility_name"))),
    _gc(getParam<Real>("gc")),
    _gcprop(coupledValue("gcprop")),
    _visco(getParam<Real>("visco")),
    _kdamage(getParam<Real>("kdamage")),

    _Hist(declareProperty<Real>("Hist")),
    _Hist_old(getMaterialPropertyOld<Real>("Hist")),
    _elastic_energy(declareProperty<Real>(getParam<MaterialPropertyName>("elastic_energy_name"))),
    _delastic_energydc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("elastic_energy_name"), coupledName("c", 0))),
    _d2elastic_energyd2c(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("elastic_energy_name"), coupledName("c", 0), coupledName("c", 0))),
    _dstress_dc(declarePropertyDerivative<Real>(_base_name + "stress", coupledName("c", 0))),
    _sigma(declareProperty<RankTwoTensor>("sigma")),
    _sigma_pressure(declareProperty<Real>("sigma_pressure")),
    _sigma_dev(declareProperty<RankTwoTensor>("sigma_dev")),
    _sigma_pos(declareProperty<RankTwoTensor>("sigma_pos")),
    _sigma_neg(declareProperty<RankTwoTensor>("sigma_neg")),
    _W(declareProperty<Real>("W")),
    _Wpos(declareProperty<Real>("Wpos")),
    _Wneg(declareProperty<Real>("Wneg")),
    //invariants for debugging
    _I1_pos(declareProperty<Real>("I1_pos")),
    _I3_pos(declareProperty<Real>("I3_pos")),
    _I1_neg(declareProperty<Real>("I1_neg")),
    _I3_neg(declareProperty<Real>("I3_neg")),
    _elastic_energy_total(declareProperty<Real>("elastic_energy_total"))
{
}

void
ArtVisJ2StressLIPIT::initialSetup()
{
  _flow_stress_material = &getMaterial("flow_stress_material");

  // Enforce isotropic elastic tensor
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("ArtVisJ2StressLIPIT requires an isotropic elasticity tensor");
}

void
ArtVisJ2StressLIPIT::initQpStatefulProperties()
{
  ComputeLagrangianStressPK1::initQpStatefulProperties();
  _be[_qp].setToIdentity();
  _ep[_qp] = 0;
  _Hist[_qp] = 0.;
}

void
ArtVisJ2StressLIPIT::computeQpPK1Stress()
{
  usingTensorIndices(i, j, k, l, m);
  const Real E = ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  const Real K = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  const auto I = RankTwoTensor::Identity();
  const auto Fit = _F[_qp].inverse().transpose();
  const auto detJ = _F[_qp].det();

  // Update configuration
  RankTwoTensor f = _inv_df[_qp].inverse();
  RankTwoTensor f_bar = f / std::cbrt(f.det());

  // Elastic predictor
  _be[_qp] = f_bar * _be_old[_qp] * f_bar.transpose();
  RankTwoTensor s = G * _be[_qp].deviatoric();
  _Np[_qp] = MooseUtils::absoluteFuzzyEqual(s.norm(), 0) ? std::sqrt(1. / 2.) * I
                                                         : std::sqrt(3. / 2.) * s / s.norm();
  Real s_eff = s.doubleContraction(_Np[_qp]);

  // Compute the derivative of the strain before return mapping
  if (_fe_problem.currentlyComputingJacobian())
    _d_be_d_F = _F_old[_qp].inverse().times<l, m, i, j, k, m>(
        (I.times<i, k, j, l>(f_bar * _be_old[_qp].transpose()) +
         I.times<j, k, i, l>(f_bar * _be_old[_qp])) /
            std::cbrt(f.det()) -
        2. / 3. * _be[_qp].times<i, j, l, k>(_inv_df[_qp]));

  // Check for plastic loading and do return mapping
  Real delta_ep = 0;
  if (computeResidual(s_eff, 0) > 0)
  {
    // Initialize the derivative of the internal variable
    if (_fe_problem.currentlyComputingJacobian())
    {
      _d_deltaep_d_betr.zero();
      if (MooseUtils::absoluteFuzzyEqual(s.norm(), 0))
        _d_n_d_be.zero();
      else
        _d_n_d_be = G / std::sqrt(6) / s.norm() *
                    (3 * I.times<i, k, j, l>(I) - 2 * _Np[_qp].times<i, j, k, l>(_Np[_qp]) -
                     I.times<i, j, k, l>(I));
    }

    returnMappingSolve(s_eff, delta_ep, _console);

    // Correct the derivative of the strain after return mapping
    if (_fe_problem.currentlyComputingJacobian())
      _d_be_d_F -=
          2. / 3. *
          (_be[_qp].trace() * _Np[_qp].times<i, j, k, l>(_d_deltaep_d_betr) +
           delta_ep * _Np[_qp].times<i, j, k, l>(I) + delta_ep * _be[_qp].trace() * _d_n_d_be) *
          _d_be_d_F;
  }

  // Update intermediate and current configurations
  _ep[_qp] = _ep_old[_qp] + delta_ep;
  _be[_qp] -= 2. / 3. * delta_ep * _be[_qp].trace() * _Np[_qp];

  //TESTING: declare plastic strain from analytic derivation
  _plastic_strain[_qp] = std::pow(detJ, - 1. / 3.) * _be[_qp].inverse() * _deformation_gradient[_qp];
  _pnorm[_qp] = std::sqrt(2. / 3.) * std::sqrt(_plastic_strain[_qp].doubleContraction(_plastic_strain[_qp]));
  _plastic_strain_rate[_qp] = (_plastic_strain[_qp] - _plastic_strain_old[_qp]) / _dt;
  _hsp[_qp] = _cauchy_stress[_qp].doubleContraction(_plastic_strain_rate[_qp]);


////////////////////////////////////////////////////////////////////////////////////

  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);

  Real I1;
  Real I3;
  RankTwoTensor B = _F[_qp] * _F[_qp].transpose(); //check
  RankTwoTensor C = B.transpose(); //check
  I1 = C.trace();
  I3 = C.det();
  Real dW_dI1 = mu / 2.;
  Real dW_dI3 = (lambda * log(std::pow(I3, 0.5)) - mu) / (2. * I3);

  //After this point, fracture is applyed as the plastic update has already been passed to the stress tensors
  ///////////////////////////FRACTURE/////////////////////////////
  //assign L and kappa
  _kappa[_qp] = _gcprop[_qp] * _l;
  _L[_qp] = 1. / (_gcprop[_qp] * _visco);

  //assign degradation derivatives
  Real factor = (_c[_qp] > 1.) ? 0. : 1.;
  Real D = factor * std::pow((1. - _c[_qp]), 2.) * (1. - _kdamage) + _kdamage;
  Real dDdc = factor * -2. * (1. - _c[_qp]) * (1. - _kdamage);
  Real d2Dd2c = factor * 2. * (1. - _kdamage);

  Real elastic_energy = (lambda / 2.) * std::pow(log(std::pow(I3, 0.5)), 2.) + (mu / 2.) * (I1 - 3.) - mu * log(std::pow(I3, 0.5));
  _elastic_energy_total[_qp] = elastic_energy; //undifferentiated elastic energy directly from cauchy tensor

  //get the eigen decompostion of the green lagrange strain tensor

  RankTwoTensor eigvec_strain;
  std::vector<Real> eigval_strain(LIBMESH_DIM);
  RankFourTensor ProjectionTensor = C.positiveProjectionEigenDecomposition(eigval_strain, eigvec_strain); //create positive projection tensor

  RankTwoTensor Cpos = ProjectionTensor * C;
  RankTwoTensor Bpos = Cpos.transpose();
  RankTwoTensor Cneg = C - Cpos;
  RankTwoTensor Bneg = Cneg.transpose();

  //get the invariants of the positive eigenvalues
  Real eigval_strain_sum = 0.;
  Real eigval_strain_mult_pos = 1.;
  Real eigval_strain_mult_neg = 1.;

  for (const auto i : make_range(Moose::dim)){
    eigval_strain_sum += eigval_strain[i];
    if (eigval_strain[i] > 0.){
      eigval_strain_mult_pos *= eigval_strain[i];
    }
    else if (eigval_strain[i] < 0.){
      eigval_strain_mult_neg *= eigval_strain[i];
    }
  }

  //get positive and negative sums
  const Real eigval_strain_sum_pos = 0.5 * (std::abs(eigval_strain_sum) + eigval_strain_sum);
  const Real eigval_strain_sum_neg = - 0.5 * (std::abs(eigval_strain_sum) - eigval_strain_sum);


  Real I1_pos = eigval_strain_sum_pos; //sum of positive eigenvalues
  _I1_pos[_qp] = I1_pos;
  Real I3_pos = eigval_strain_mult_pos; //product of positive eigenvalues
  _I3_pos[_qp] = I3_pos;

  Real I1_neg = eigval_strain_sum_neg; //sum of negative eigenvalues
  _I1_neg[_qp] = I1_neg;
  Real I3_neg = eigval_strain_mult_neg; //product of negative eigenvalues
  _I3_neg[_qp] = I3_neg;

  Real elastic_energy_pos = (lambda / 2.) * std::pow(log(std::pow(I3_pos, 0.5)), 2.) + (mu / 2.) * (I1_pos - 3.) - mu * log(std::pow(I3_pos, 0.5));
  Real elastic_energy_neg = elastic_energy - elastic_energy_pos;

  //define total elastic energy and history variable

  _W[_qp] = D * elastic_energy_pos + elastic_energy_neg;

  if (elastic_energy_pos > _Hist_old[_qp]){
    _Hist[_qp] = elastic_energy_pos;
  }else{
    _Hist[_qp] = _Hist_old[_qp];
  }

  _Wpos[_qp] = elastic_energy_pos;
  _Wneg[_qp] = elastic_energy_neg;

  //penalize the stress with the degradation value, currently penalizing the entire stress tensor
  //penalize by vol-dev decomposition

  //_pk1_stress[_qp] *= D; //JUST CHANGED; CHECK AND FIX IF NECESSARY
  _dstress_dc[_qp] = D;
  _elastic_energy[_qp] = (D * _Hist[_qp]) + (_gcprop[_qp] * std::pow(_c[_qp], 2.) / (2. * _l)); //sum of penalized elastic energy plus fractured new surface energy
  _delastic_energydc[_qp] = (dDdc * _Hist[_qp]) + (_gcprop[_qp] * _c[_qp] / _l);
  _d2elastic_energyd2c[_qp] = (d2Dd2c * _Hist[_qp]) + (_gcprop[_qp] / _l);


/////////////////////////////////////////////////////////////////////////////

  s = G * _be[_qp].deviatoric();
  //RankTwoTensor tau = (K / 2.) * (std::pow(_deformation_gradient[_qp].det(), 2.) - 1.) * I + s;
  //TESTING: definition of strain energy and derivatives
  
  RankTwoTensor sigma = (2. / std::pow(I3, 0.5)) * ((I3 * dW_dI3 * I) + (dW_dI1 * B));
  
  //get positive and negative stress tensors
  //NOTE: we use the spectral decomposition of the green lagrange strain tensor to get the stress too

  RankTwoTensor sigma_pos = ProjectionTensor * sigma; //positive part of stress tensor
  RankTwoTensor sigma_neg = sigma - sigma_pos; //negative part of stress tensor

  //output stress tensors
  RankTwoTensor sigma_penalized = D * sigma_pos + sigma_neg;
  RankTwoTensor tau = detJ * sigma_penalized;
  _sigma[_qp] = sigma_penalized;
  _sigma_pressure[_qp] = - (1. / 3.) * sigma_penalized.trace();
  _sigma_dev[_qp] = sigma_penalized.deviatoric();

  _sigma_pos[_qp] = sigma_pos;
  _sigma_neg[_qp] = sigma_neg;
  

  _pk1_stress[_qp] = detJ * (sigma_penalized) * Fit;

  /////////////////////////ARTIFICIAL VISCOSITY////////////////////////////
  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour

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
  _pk1_stress[_qp] += P_av * I2;

  // Compute the consistent tangent, i.e. the derivative of the PK1 stress w.r.t. the deformation
  // gradient.
  if (_fe_problem.currentlyComputingJacobian())
  {
    RankFourTensor d_tau_d_F = K * detJ * detJ * I.times<i, j, k, l>(Fit) +
                               G * (_d_be_d_F - I.times<i, j, k, l>(I) * _d_be_d_F / 3);
    _pk1_jacobian[_qp] = Fit.times<m, j, i, m, k, l>(d_tau_d_F) - Fit.times<k, j, i, l>(tau * Fit);
  }

}

Real
ArtVisJ2StressLIPIT::computeReferenceResidual(const Real & effective_trial_stress,
                                                              const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  return effective_trial_stress - G * scalar * _be[_qp].trace();
}

Real
ArtVisJ2StressLIPIT::computeResidual(const Real & effective_trial_stress,
                                                     const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return effective_trial_stress - G * scalar * _be[_qp].trace() - _H[_qp];
}

Real
ArtVisJ2StressLIPIT::computeDerivative(const Real & /*effective_trial_stress*/,
                                                       const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return -G * _be[_qp].trace() - _dH[_qp];
}

void
ArtVisJ2StressLIPIT::preStep(const Real & scalar, const Real & R, const Real & J)
{
  if (!_fe_problem.currentlyComputingJacobian())
    return;

  const auto I = RankTwoTensor::Identity();
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  _d_R_d_betr =
      G * _Np[_qp] - G * scalar * I - (G * _be[_qp].trace() + _dH[_qp]) * _d_deltaep_d_betr;
  _d_J_d_betr = -G * I - _d2H[_qp] * _d_deltaep_d_betr;
  _d_deltaep_d_betr += -1 / J * _d_R_d_betr + R / J / J * _d_J_d_betr;

}