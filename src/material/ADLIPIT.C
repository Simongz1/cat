#include "ADLIPIT.h"

registerMooseObject("mlApp", ADLIPIT);

InputParameters
ADLIPIT::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADComputeStressBase>::validParams();
  params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("Finite fracture - plasticity");
  params.addParam<MaterialPropertyName>("elasticity_tensor", "elasticity_tensor", "The name of the elasticity tensor.");
  params.addRequiredParam<MaterialName>("flow_stress_material", "The material defining the flow stress");
  /////////////
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredParam<Real>("element_size", "element_size");

  //required variables
  params.addRequiredCoupledVar("c", "fracture variable");
  params.addRequiredCoupledVar("gc", "gc");
  params.addRequiredCoupledVar("displacements", "displacements");

  params.addParam<bool>("use_custom", true, "use custom deformation gradient calculation");
  return params;
}

ADLIPIT::ADLIPIT(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ADComputeStressBase>(parameters),
    ADSingleVariableReturnMappingSolution(parameters),
    _elasticity_tensor_name(getParam<MaterialPropertyName>("elasticity_tensor")),
    _elasticity_tensor(getADMaterialProperty<RankFourTensor>(_elasticity_tensor_name)),
    _F(declareADProperty<RankTwoTensor>("deformation_gradient")),
    _F_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _Fhat(declareADProperty<RankTwoTensor>("Fhat")),

    _ep_name("ep"),
    _ep(declareADProperty<Real>(_ep_name)),
    _ep_old(getMaterialPropertyOldByName<Real>(_ep_name)),
    _ep_dot(declareADProperty<Real>("ep_dot")),
    _be(declareADProperty<RankTwoTensor>("volume_preserving_elastic_left_cauchy_green_strain")),
    _be_old(getMaterialPropertyOldByName<RankTwoTensor>("volume_preserving_elastic_left_cauchy_green_strain")),
    _Np(declareADProperty<RankTwoTensor>("flow_direction")),

    //treating Fp as a stateful property
    _Fp(declareADProperty<RankTwoTensor>("Fp")),
    _Fp_old(getMaterialPropertyOld<RankTwoTensor>("Fp")),
    _Fe(declareADProperty<RankTwoTensor>("Fe")),
    _Fe_old(getMaterialPropertyOld<RankTwoTensor>("Fe")),

    //generate strains and rates

    _Ee(declareADProperty<RankTwoTensor>("Ee")),
    _Ee_dot(declareADProperty<RankTwoTensor>("Ee_dot")),

    _Ep(declareADProperty<RankTwoTensor>("Ep")),
    _Ep_dot(declareADProperty<RankTwoTensor>("Ep_dot")),

    _flow_stress_material(nullptr),
    _flow_stress_name("flow_stress"),

    _H(getADMaterialPropertyByName<Real>(_flow_stress_name)),
    _dH(getMaterialPropertyDerivativeByName<Real>(_flow_stress_name, _ep_name)),
    _d2H(getMaterialPropertyDerivativeByName<Real>(_flow_stress_name, _ep_name, _ep_name)),

    /////
    _rho(getADMaterialProperty<Real>("density")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("element_size")),

    /////////////////

    //request fracture stuff
    _c(adCoupledValue("c")),
    _gc(adCoupledValue("gc")),

    _Hist(declareADProperty<Real>("Hist")),
    _Hist_old(getMaterialPropertyOld<Real>("Hist")),
  
    _W0(declareADProperty<Real>("W0")),
    _W(declareADProperty<Real>("W")),
    _Wpos(declareADProperty<Real>("Wpos")),
    _Wneg(declareADProperty<Real>("Wneg")),

    //invariants for debugging
    _inv_Cp(declareADProperty<RankTwoTensor>("inv_Cp")),

    _Cp(declareADProperty<RankTwoTensor>("Cp")),
    _inv_Cp_old(getMaterialPropertyOld<RankTwoTensor>("inv_Cp")),

    _Cp_old(getMaterialPropertyOld<RankTwoTensor>("Cp")),
    _Ce(declareADProperty<RankTwoTensor>("Ce")),

    _Ce_old(getMaterialPropertyOld<RankTwoTensor>("Ce")),
    _Ep_old(getMaterialPropertyOld<RankTwoTensor>("Ep")),
    _Ee_old(getMaterialPropertyOld<RankTwoTensor>("Ee")),

    _S(declareADProperty<RankTwoTensor>("S")),
    _HS_elastic(declareADProperty<Real>("HS_elastic")),
    _HS_plastic(declareADProperty<Real>("HS_plastic")),

    _D(declareADProperty<Real>("D")),
    _kdamage(getADMaterialProperty<Real>("kdamage")),
    _Je(declareADProperty<Real>("Je")),
    _Jp(declareADProperty<Real>("Jp")),
    _J(declareADProperty<Real>("J")),
    _Fres(declareADProperty<RankTwoTensor>("Fres")),
    _E(declareADProperty<RankTwoTensor>("E")),
    _E_dot(declareADProperty<RankTwoTensor>("E_dot")),

    //cauchy stress
    _sigma(declareADProperty<RankTwoTensor>("sigma")),
    _ndisp(coupledComponents("displacements")),

    //get base increments
    _strain_increment(getADMaterialProperty<RankTwoTensor>("strain_increment")),
    _rotation_increment(getADMaterialProperty<RankTwoTensor>("rotation_increment")),
    _F_from_increment(declareADProperty<RankTwoTensor>("F_from_increment")),
    _use_custom(getParam<bool>("use_custom"))
{
  //extract displacements
  _grad_disp.reserve(_ndisp);
  _grad_disp_old.reserve(_ndisp);

  for (unsigned int i = 0; i < _ndisp; ++i){
    _grad_disp.push_back(&adCoupledGradient("displacements", i));
    _grad_disp_old.push_back(&coupledGradientOld("displacements", i));
  }
}

void
ADLIPIT::initialSetup()
{
  _flow_stress_material = &getMaterial("flow_stress_material");
}

void
ADLIPIT::initQpStatefulProperties()
{
  ADComputeStressBase::initQpStatefulProperties();
  _be[_qp].setToIdentity();
  _ep[_qp] = 0;
  _Fp[_qp].setToIdentity();
  _Fe[_qp].setToIdentity();
  _Cp[_qp].setToIdentity();
  _F[_qp].setToIdentity();
  _Hist[_qp] = 0;
  _Fhat[_qp].setToIdentity();
}

void
ADLIPIT::computeQpStress()
{
  ADRankTwoTensor I2;
  I2.setToIdentity();
  const ADReal G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  const ADReal K = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  const ADRankTwoTensor I = ADRankTwoTensor::Identity();

  //compute AD version of the incremental deformation gradient
  //form a tensor with rows occupied by displacements

  std::vector<ADRealGradient> grad_disp_vect(_ndisp);

  for (unsigned int i = 0; i < _ndisp; i++){
    grad_disp_vect[i] = (*_grad_disp[i])[_qp];
  }

  ADRankTwoTensor A = ADRankTwoTensor::initializeFromRows(
    grad_disp_vect[0], grad_disp_vect[1], grad_disp_vect[2]);

  //old deformation gradient
  const auto & g0_old = (*_grad_disp_old[0])[_qp];
  const auto & g1_old = (*_grad_disp_old[1])[_qp];
  const auto & g2_old = (*_grad_disp_old[2])[_qp];

  ADRealGradient g0_old_ad(g0_old(0), g0_old(1), g0_old(2));
  ADRealGradient g1_old_ad(g1_old(0), g1_old(1), g1_old(2));
  ADRealGradient g2_old_ad(g2_old(0), g2_old(1), g2_old(2));

  ADRankTwoTensor Fbar = ADRankTwoTensor::initializeFromRows(
    g0_old_ad, g1_old_ad, g2_old_ad);

  //update A (incremental deformation gradient)
  A -= Fbar;

  //update Fbar (old deformation gradient)
  Fbar.addIa(1.);

  ADRankTwoTensor F_incremental;
  ADRankTwoTensor Fit;
  ///FORM CASES FOR UPDATE
  if (_use_custom){
    //custom branch
    F_incremental = A * Fbar.inverse();
    F_incremental.addIa(1.);

  }else{
    //from ADComputeFiniteStrain
    F_incremental = _rotation_increment[_qp] * (_strain_increment[_qp] + I);
    _F_from_increment[_qp] = F_incremental;
  }

  //after deciding method, compute other stuff
  _F[_qp] = F_incremental * Fbar;
  _Fhat[_qp] = F_incremental;
  Fit = _F[_qp].inverse().transpose();

  ////////////////////////////////////
  ADRankTwoTensor f = _Fhat[_qp];
  ADReal Jhat = f.det();
  ADRankTwoTensor f_bar = f / MetaPhysicL::cbrt(Jhat);
  ADReal J = _F[_qp].det();

  // Elastic predictor
  _be[_qp] = f_bar * _be_old[_qp] * f_bar.transpose();
  ADRankTwoTensor s = G * _be[_qp].deviatoric();
  ADReal snorm = MetaPhysicL::sqrt(s.doubleContraction(s));
  _Np[_qp] = MooseUtils::absoluteFuzzyEqual(snorm, ADReal(0)) ? std::sqrt(1. / 2.) * I
                                                         : std::sqrt(3. / 2.) * s / snorm;
  ADReal s_eff = s.doubleContraction(_Np[_qp]);

  // Check for plastic loading and do return mapping
  ADReal delta_ep = 0;
  if (MetaPhysicL::raw_value(computeResidual(s_eff, 0)) > 0)
  {
    returnMappingSolve(s_eff, delta_ep, _console);
  }

  // Update intermediate and current configurations
  _ep[_qp] = _ep_old[_qp] + delta_ep;
  _be[_qp] -= 2. / 3. * delta_ep * _be[_qp].trace() * _Np[_qp];

  //obtain inverse plastic volume preserving C tensor

  //compute F_bar at n+1
  //RankTwoTensor F = f * _F_old[_qp].inverse(); //equivalent to computing F[n+1] = f[n+1]F[n].inverse();
  ADRankTwoTensor F = _F[_qp];
  ADReal detJ = F.det();
  ADRankTwoTensor F_bar = MetaPhysicL::pow(J, - 1. / 3.) * F;

  //use identity to get volume preserving C^p^-1

  _inv_Cp[_qp] = F_bar.inverse() * _be[_qp] * F_bar.inverse().transpose();
  _Cp[_qp] = _inv_Cp[_qp].inverse();

  ADRankTwoTensor V;
  std::vector<ADReal> lam(3);
  
  ADRankTwoTensor Csym = 0.5 * (_Cp[_qp] + _Cp[_qp].transpose());
  Csym.symmetricEigenvaluesEigenvectors(lam, V);

  ADRankTwoTensor diag;
  diag.zero();
  for (unsigned int i = 0; i < 3; ++i){
    diag(i,i) = MetaPhysicL::sqrt(std::max(lam[i], ADReal(1e-10)));
  }

  //compute Fp
  _Fp[_qp] = V * diag * V.transpose();

  //compute elastic deformation gradient
  _Fe[_qp] = _F[_qp] * _Fp[_qp].inverse();
  _Ce[_qp] = _Fe[_qp].transpose() * _Fe[_qp];

  //compute eleastic green strain
  _Ee[_qp] = 0.5 * (_Fe[_qp].transpose() * _Fe[_qp] - I2);

  //use this to compute plastic strain
  _Ep[_qp] = 0.5 * (_Cp[_qp] - I);

  //compute deformation gradient rates and strain rates
  ADRankTwoTensor F_dot, Fe_dot, Fp_dot;
  Fp_dot = (1. / _dt) * (_Fp[_qp] - _Fp_old[_qp]);
  Fe_dot = (1. / _dt) * (_Fe[_qp] - _Fe_old[_qp]);
  F_dot = (1. / _dt) * (_F[_qp] - _F_old[_qp]);

  _Ep_dot[_qp] = 0.5 * (Fp_dot.transpose() * _Fp[_qp] + _Fp[_qp].transpose() * Fp_dot);
  _Ee_dot[_qp] = 0.5 * (Fe_dot.transpose() * _Fe[_qp] + _Fe[_qp].transpose() * Fe_dot);
  _ep_dot[_qp] = (1. / _dt) * (_ep[_qp] - _ep_old[_qp]);

  ///invariants for elastic energy calculation

  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);

  //elastic energy stuff for fracture

  ////////////////////////////////////
  //compute the positive and negative invariants for strain energy
  auto braket_pos = [](ADReal val) {return val > 0. ? val : 0.;};
  auto braket_neg = [](ADReal val) {return val < 0. ? val : 0.;};

  //cache stuff here
  ADReal Je = _Fe[_qp].det();
  ADReal trbe = _be[_qp].trace();
  ADReal x = trbe - 3.;

  //define positive and negative energies based on their respective variables
  ADReal Wvolpos, Wvolneg, Wdevpos, Wdevneg;
  ADReal Jpos, Jneg, xpos, xneg;

  //compute indicators
  Jpos = 1. + braket_pos(J - 1.);
  Jneg = 1. + braket_neg(J - 1.);
  xpos = braket_pos(trbe - 3.);
  xneg = braket_neg(trbe - 3.);

  //compute energies
  Wvolpos = computeWvol(Jpos); Wvolneg = computeWvol(Jneg);
  Wdevpos = computeWinv(xpos); Wdevneg = computeWinv(xneg);

  //compute degradation
  _D[_qp] = (1 - _kdamage[_qp]) * MetaPhysicL::pow(1 - _c[_qp], 2. * (1. + _ep[_qp])) + _kdamage[_qp];

  //compose penalized strain energy
  _W[_qp] = (_D[_qp] * Wvolpos) + Wvolneg + (_D[_qp] * Wdevpos) + Wdevneg; 
  _Wpos[_qp] = Wvolpos + Wdevpos;
  _Wneg[_qp] = Wvolneg + Wdevneg;

  //compute history variable

  if (_Wpos[_qp] > _Hist_old[_qp]){
    _Hist[_qp] = _Wpos[_qp];
  }else{
    _Hist[_qp] = _Hist_old[_qp];
  }

  //compute split kirchhoff stress: volumetric part
  ADReal dWdJpos, dJposdJ;
  ADReal dWdJneg, dJnegdJ;

  dWdJpos = 0.5 * K * (Jpos - (1. / Jpos)); dJposdJ = (Je > 1) ? 1. : 0.;
  dWdJneg = 0.5 * K * (Jneg - (1. / Jneg)); dJnegdJ = (Je < 1) ? 1. : 0.;

  //define volumetric kirchhoff
  ADRankTwoTensor tauvol;
  tauvol = _F[_qp].det() * (_D[_qp] * dWdJpos * dJposdJ + dWdJneg * dJnegdJ) * I2;

  //define deviatoric tau and degrade
  ADRankTwoTensor taudev;
  taudev = _D[_qp] * mu * _be[_qp].deviatoric();

  ADRankTwoTensor tau = tauvol + taudev;

  ADRankTwoTensor pk1 = tau * _F[_qp].inverse().transpose();
  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour

  ADReal ss = MetaPhysicL::sqrt(K / _rho[_qp]);
	
  //Compute artificial viscosity term
  ADReal P_av;
  ADReal Je_dot;
  Je_dot = ((_Fe[_qp].det() - _Fe_old[_qp].det()) / _dt);
  ADReal Jtot = _F[_qp].det();

  P_av = _C0 * _rho[_qp] * (Je_dot * MetaPhysicL::abs(Je_dot) / MetaPhysicL::pow(Je, 2.0)) * std::pow(_Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (Je_dot / Je) * _Le;
  //pk1 += P_av * I;

  //compute here PK2 stress
  _S[_qp] = _F[_qp].inverse() * pk1 + P_av * I;

  //cauchy
  _stress[_qp] = (1. / J) * tau + P_av * I;
  _Je[_qp] = _Fe[_qp].det();
  _Jp[_qp] = _Fp[_qp].det();
  _J[_qp] = J;

  _Fres[_qp] = _F[_qp] - _Fe[_qp] * _Fp[_qp];
  _E_dot[_qp] = 0.5 * (F_dot.transpose() * F + F.transpose() * F_dot);
  _E[_qp] = 0.5 * (_F[_qp].transpose() * _F[_qp] - I2);

  // Compute the consistent tangent, i.e. the derivative of the PK1 stress w.r.t. the deformation
  // gradient.

  _HS_elastic[_qp] = P_av * _Ee_dot[_qp].trace();
  _HS_plastic[_qp] = _S[_qp].doubleContraction(_Ep_dot[_qp]);
}

Real
ADLIPIT::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                              const ADReal & scalar)
{
  const ADReal G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  return MetaPhysicL::raw_value(effective_trial_stress - G * scalar * _be[_qp].trace());
}

ADReal
ADLIPIT::computeResidual(const ADReal & effective_trial_stress,
                                                     const ADReal & scalar)
{
  const ADReal G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return (effective_trial_stress - G * scalar * _be[_qp].trace() - _H[_qp]);
}

ADReal
ADLIPIT::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                       const ADReal & scalar)
{
  const ADReal G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return (-G * _be[_qp].trace() - _dH[_qp]);
}

//define a function to get strain energy
ADReal
ADLIPIT::computeWvol(const ADReal & J)
{
  const ADReal K = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);

  //evaluate expression
  //Real Winv = lambda * ((std::pow(I3, 2.) - 1.) / 4.) - ((lambda / 2.) - mu) * std::log(I3);
  ADReal Wvol = (K / 4.) * (J * J - 1. - 2. * MetaPhysicL::log(J));
  return Wvol;
}

ADReal
ADLIPIT::computeWinv(const ADReal & x)
{
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  ADReal Winv = 0.5 * mu * x;
  return Winv;
}