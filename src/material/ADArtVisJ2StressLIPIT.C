//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADArtVisJ2StressLIPIT.h"

registerMooseObject("catApp", ADArtVisJ2StressLIPIT);

InputParameters
ADArtVisJ2StressLIPIT::validParams()
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
  params.addRequiredParam<Real>("ep_ref", "refference effective plastic strain");
  return params;
}

ADArtVisJ2StressLIPIT::ADArtVisJ2StressLIPIT(
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
    //treating Fp as a stateful property
    _Fp(declareProperty<RankTwoTensor>("Fp")),
    _Fp_old(getMaterialPropertyOld<RankTwoTensor>("Fp")),
    _Fe(declareProperty<RankTwoTensor>("Fe")),
    _Fe_old(getMaterialPropertyOld<RankTwoTensor>("Fe")),

    //generate strains and rates

    _Ee(declareProperty<RankTwoTensor>("Ee")),
    _Ee_dot(declareProperty<RankTwoTensor>("Ee_dot")),

    _Ep(declareProperty<RankTwoTensor>("Ep")),
    _Ep_dot(declareProperty<RankTwoTensor>("Ep_dot")),

    _flow_stress_material(nullptr),
    _flow_stress_name(_base_name + "flow_stress"),
    _H(getMaterialPropertyByName<Real>(_flow_stress_name)),
    _dH(getMaterialProperty<Real>("dH")),
    _d2H(getMaterialProperty<Real>("d2H")),
    /////
    _rho(getADMaterialProperty<Real>("density")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("element_size")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _cauchy_stress(getMaterialProperty<RankTwoTensor>("cauchy_stress")),

    /////////////////

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
    _elastic_energy_total(declareProperty<Real>("elastic_energy_total")),
    _ep_ref(getParam<Real>("ep_ref"))
{
}

void
ADArtVisJ2StressLIPIT::initialSetup()
{
  _flow_stress_material = &getMaterial("flow_stress_material");

  // Enforce isotropic elastic tensor
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("ADArtVisJ2StressLIPIT requires an isotropic elasticity tensor");
}

void
ADArtVisJ2StressLIPIT::initQpStatefulProperties()
{
  ComputeLagrangianStressPK1::initQpStatefulProperties();
  _be[_qp].setToIdentity();
  _ep[_qp] = 0;
  _Fp[_qp].setToIdentity(); //stateful plastic deformation gradient
  _Fe[_qp].setToIdentity();
}

void
ADArtVisJ2StressLIPIT::computeQpPK1Stress()
{
  usingTensorIndices(i, j, k, l, m);
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

  //incrementally update Fp and Fe
  RankTwoTensor delta_Fp = I + delta_ep * _Np[_qp]; //check
  _Fp[_qp] = delta_Fp * _Fp_old[_qp]; //check
  _Fe[_qp] = _F[_qp] * _Fp[_qp].inverse(); //check

  //generate rates for visualization
  _Ee[_qp] = (1. / 2.) * ((_Fe[_qp].transpose() * _Fe[_qp]) - I);
  RankTwoTensor Fe_dot = (_Fe[_qp] - _Fe_old[_qp]) / _dt;
  _Ee_dot[_qp] = (1. / 2.) * ((Fe_dot.transpose() * _Fe[_qp]) + _Fe[_qp].transpose() * Fe_dot);

  _Ep[_qp] = (1. / 2.) * ((_Fp[_qp].transpose() * _Fp[_qp]) - I);
  RankTwoTensor Fp_dot = (_Fp[_qp] - _Fp_old[_qp]) / _dt;
  _Ep_dot[_qp] = (1. / 2.) * ((Fp_dot.transpose() * _Fp[_qp]) + _Fp[_qp].transpose() * Fp_dot);


  ///invariants for elastic energy calculation

  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);

  Real I1;
  Real I3;
  RankTwoTensor B = _Fe[_qp] * _Fe[_qp].transpose();
  RankTwoTensor C = B.transpose(); //check
  I1 = C.trace();
  I3 = C.det();
  Real dW_dI1 = mu / 2.;
  Real dW_dI3 = (lambda * log(std::pow(I3, 0.5)) - mu) / (2. * I3);

  //elastic energy stuff for fracture

  _kappa[_qp] = _gcprop[_qp] * _l;
  _L[_qp] = 1. / (_gcprop[_qp] * _visco);

  //assign degradation derivatives
  Real S = 1. - _c[_qp]; //degradation value
  Real e_norm = _ep[_qp] / _ep_ref; //normalized plastic strain
  Real D = std::pow(S, 2 * e_norm) + _kdamage;
  Real dDdc = - 2. * e_norm * std::pow(S, (2.* e_norm) - 1.); //derivative of degradation w.r.t. c
  Real d2Dd2c = 2. * e_norm * (2. * e_norm - 1.) * std::pow(S, (2. * e_norm) - 2.);

  Real elastic_energy = (lambda / 2.) * std::pow(log(std::pow(I3, 0.5)), 2.) + (mu / 2.) * (I1 - 3.) - mu * log(std::pow(I3, 0.5));
  _elastic_energy_total[_qp] = elastic_energy; //undifferentiated elastic energy directly from cauchy tensor


  //stress stuff
  RankTwoTensor sigma = (2. / std::pow(I3, 0.5)) * ((I3 * dW_dI3 * I) + (dW_dI1 * B)); //total stress tensor

  //get the strain energy where volumetric strain is positive

  Real W_pos;
  Real W_neg;
  RankTwoTensor sigma_pos;
  RankTwoTensor sigma_neg;

  //get positive loading regions, ORIGINALLY JUST WHERE J > 1

  //write stress positive and negative
  if (_F[_qp].det() > 1.){
    W_pos = elastic_energy;
    W_neg = 0;
    sigma_pos = sigma;
    sigma_neg = sigma - sigma_pos;
    _sigma_pos[_qp] = sigma_pos;
    _sigma_neg[_qp] = sigma_neg;
    _Wpos[_qp] = W_pos;
    _Wneg[_qp] = W_neg;
  }

  if(_F[_qp].det() <= 1.){
    W_neg = elastic_energy;
    W_pos = 0;
    sigma_neg = sigma;
    sigma_pos = sigma - sigma_neg;
    _sigma_pos[_qp] = sigma_pos;
    _sigma_neg[_qp] = sigma_neg;
    _Wpos[_qp] = W_pos;
    _Wneg[_qp] = W_neg;
  }

  //define total elastic energy and history variable

  _W[_qp] = D * _Wpos[_qp] + _Wneg[_qp];

  if (_Wpos[_qp] > _Hist_old[_qp]){
    _Hist[_qp] = _Wpos[_qp];
  }else{
    _Hist[_qp] = _Hist_old[_qp];
  }

  //damage stuff

  _dstress_dc[_qp] = D;
  _elastic_energy[_qp] = (D * _Hist[_qp]) + (_gcprop[_qp] * std::pow(_c[_qp], 2.) / (2. * _l)); //sum of penalized elastic energy plus fractured new surface energy
  _delastic_energydc[_qp] = (dDdc * _Hist[_qp]) + (_gcprop[_qp] * _c[_qp] / _l);
  _d2elastic_energyd2c[_qp] = (d2Dd2c * _Hist[_qp]) + (_gcprop[_qp] / _l);
  
  //output stress tensors
  RankTwoTensor sigma_penalized = D * sigma_pos + sigma_neg;
  RankTwoTensor tau = _F[_qp].det() * sigma_penalized;
  _sigma[_qp] = sigma_penalized;
  _sigma_pressure[_qp] = - (1. / 3.) * sigma_penalized.trace();
  _sigma_dev[_qp] = sigma_penalized.deviatoric();

  _pk1_stress[_qp] = detJ * (sigma_penalized) * Fit;

  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour

  Real ss = std::sqrt(K / _rho[_qp].value());
	
  //Compute artificial viscosity term
  Real P_av;
  Real Je;
  Real Je_dot;
  Je = _deformation_gradient[_qp].det();
  Je_dot = ((_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det()) / _dt);

  P_av = _C0 * _rho[_qp].value() * (Je_dot * std::abs(Je_dot) / std::pow(Je, 2.0)) * std::pow(_Le, 2.0);
  P_av += _C1 * _rho[_qp].value() * ss * (Je_dot / Je) * _Le;
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
ADArtVisJ2StressLIPIT::computeReferenceResidual(const Real & effective_trial_stress,
                                                              const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  return effective_trial_stress - G * scalar * _be[_qp].trace();
}

Real
ADArtVisJ2StressLIPIT::computeResidual(const Real & effective_trial_stress,
                                                     const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return effective_trial_stress - G * scalar * _be[_qp].trace() - _H[_qp];
}

Real
ADArtVisJ2StressLIPIT::computeDerivative(const Real & /*effective_trial_stress*/,
                                                       const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return -G * _be[_qp].trace() - _dH[_qp];
}

void
ADArtVisJ2StressLIPIT::preStep(const Real & scalar, const Real & R, const Real & J)
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