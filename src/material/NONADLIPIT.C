#include "NONADLIPIT.h"

registerMooseObject("mlApp", NONADLIPIT);

InputParameters
NONADLIPIT::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ComputeLagrangianStressPK1>::validParams();
  params += SingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("Finite fracture - plasticity");
  params.addParam<MaterialPropertyName>("elasticity_tensor", "elasticity_tensor", "The name of the elasticity tensor.");
  params.addRequiredParam<MaterialName>("flow_stress_material", "The material defining the flow stress");
  /////////////
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredParam<Real>("element_size", "element_size");

  //fracture stuff
  params.addRequiredCoupledVar("c", "fracture variable");
  params.addRequiredCoupledVar("gc", "gc");
  return params;
}

NONADLIPIT::NONADLIPIT(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ComputeLagrangianStressPK1>(parameters),
    GuaranteeConsumer(this),
    SingleVariableReturnMappingSolution(parameters),
    _elasticity_tensor_name(getParam<MaterialPropertyName>("elasticity_tensor")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>(_elasticity_tensor_name)),
    _F(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _F_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _ep_name("ep"),
    _ep(declareProperty<Real>(_ep_name)),
    _ep_old(getMaterialPropertyOldByName<Real>(_ep_name)),
    _ep_dot(declareProperty<Real>("ep_dot")),
    _be(declareProperty<RankTwoTensor>("volume_preserving_elastic_left_cauchy_green_strain")),
    _be_old(getMaterialPropertyOldByName<RankTwoTensor>("volume_preserving_elastic_left_cauchy_green_strain")),
    _Np(declareProperty<RankTwoTensor>("flow_direction")),

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
    _flow_stress_name("flow_stress"),

    _H(getMaterialPropertyByName<Real>(_flow_stress_name)),
    _dH(getDefaultMaterialPropertyByName<Real, false>(derivativePropertyName(_flow_stress_name, {_ep_name}))),
    _d2H(getDefaultMaterialPropertyByName<Real, false>(derivativePropertyName(_flow_stress_name, {_ep_name, _ep_name}))),

    /////
    _rho(getMaterialProperty<Real>("density")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("element_size")),

    /////////////////

    //request fracture stuff
    _c(coupledValue("c")),
    _gc(coupledValue("gc")),

    _Hist(declareProperty<Real>("Hist")),
    _Hist_old(getMaterialPropertyOld<Real>("Hist")),
  
    _W0(declareProperty<Real>("W0")),
    _W(declareProperty<Real>("W")),
    _Wpos(declareProperty<Real>("Wpos")),
    _Wneg(declareProperty<Real>("Wneg")),

    //invariants for debugging
    _inv_Cp(declareProperty<RankTwoTensor>("inv_Cp")),

    _Cp(declareProperty<RankTwoTensor>("Cp")),
    _inv_Cp_old(getMaterialPropertyOld<RankTwoTensor>("inv_Cp")),

    _Cp_old(getMaterialPropertyOld<RankTwoTensor>("Cp")),
    _Ce(declareProperty<RankTwoTensor>("Ce")),

    _Ce_old(getMaterialPropertyOld<RankTwoTensor>("Ce")),
    _Ep_old(getMaterialPropertyOld<RankTwoTensor>("Ep")),
    _Ee_old(getMaterialPropertyOld<RankTwoTensor>("Ee")),

    _S(declareProperty<RankTwoTensor>("S")),
    _HS_elastic(declareProperty<Real>("HS_elastic")),
    _HS_plastic(declareProperty<Real>("HS_plastic")),

    _D(declareProperty<Real>("D")),
    _pk1_pos(declareProperty<RankTwoTensor>("pk1_pos")),
    _pk1_neg(declareProperty<RankTwoTensor>("pk1_neg")),

    _kdamage(getMaterialProperty<Real>("kdamage")),
    _Je(declareProperty<Real>("Je")),
    _Jp(declareProperty<Real>("Jp")),
    _Fres(declareProperty<RankTwoTensor>("Fres")),
    _E(declareProperty<RankTwoTensor>("E")),
    _E_dot(declareProperty<RankTwoTensor>("E_dot")),

    //cauchy stress
    _sigma(declareProperty<RankTwoTensor>("sigma"))
{
}

void
NONADLIPIT::initialSetup()
{
  _flow_stress_material = &getMaterial("flow_stress_material");

  // Enforce isotropic elastic tensor
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("NONADLIPIT requires an isotropic elasticity tensor");
}

void
NONADLIPIT::initQpStatefulProperties()
{
  ComputeLagrangianStressPK1::initQpStatefulProperties();
  _be[_qp].setToIdentity();
  _ep[_qp] = 0;
  _Fp[_qp].setToIdentity();
  _Fe[_qp].setToIdentity();
  //_inv_Cp[_qp].setToIdentity();
  _Cp[_qp].setToIdentity();
  _Hist[_qp] = 0;
}

void
NONADLIPIT::computeQpPK1Stress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
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

  //obtain inverse plastic volume preserving C tensor

  //compute F_bar at n+1
  //RankTwoTensor F = f * _F_old[_qp].inverse(); //equivalent to computing F[n+1] = f[n+1]F[n].inverse();
  RankTwoTensor F = _F[_qp];
  
  RankTwoTensor F_bar = std::pow(F.det(), - 1. / 3.) * F;

  //use identity to get volume preserving C^p^-1

  _inv_Cp[_qp] = F_bar.inverse() * _be[_qp] * F_bar.inverse().transpose();
  _Cp[_qp] = _inv_Cp[_qp].inverse();

  RankTwoTensor V;
  std::vector<Real> lam(3);
  
  RankTwoTensor Csym = 0.5 * (_Cp[_qp] + _Cp[_qp].transpose());
  Csym.symmetricEigenvaluesEigenvectors(lam, V);

  RankTwoTensor diag;
  diag.zero();
  for (unsigned int i = 0; i < 3; ++i){
    diag(i,i) = std::sqrt(std::max(lam[i], 1e-10));
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
  RankTwoTensor F_dot, Fe_dot, Fp_dot;
  Fp_dot = (1. / _dt) * (_Fp[_qp] - _Fp_old[_qp]);
  Fe_dot = (1. / _dt) * (_Fe[_qp] - _Fe_old[_qp]);
  F_dot = (1. / _dt) * (_F[_qp] - _F_old[_qp]);

  _Ep_dot[_qp] = 0.5 * (Fp_dot.transpose() * _Fp[_qp] + _Fp[_qp].transpose() * Fp_dot);
  _Ee_dot[_qp] = 0.5 * (Fe_dot.transpose() * _Fe[_qp] + _Fe[_qp].transpose() * Fe_dot);
  _ep_dot[_qp] = (1. / _dt) * (_ep[_qp] - _ep_old[_qp]);

  ///invariants for elastic energy calculation

  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);

  //elastic energy stuff for fracture

  ////////////////////////////////////
  //compute the positive and negative invariants for strain energy
  auto braket_pos = [](Real val) {return val > 0. ? val : 0.;};
  auto braket_neg = [](Real val) {return val < 0. ? val : 0.;};

  //cache stuff here
  Real Je = _Fe[_qp].det();
  Real J = _F[_qp].det();
  Real trbe = _be[_qp].trace();
  Real x = trbe - 3.;

  //define positive and negative energies based on their respective variables
  Real Wvolpos, Wvolneg, Wdevpos, Wdevneg;
  Real Jpos, Jneg, xpos, xneg;

  //compute indicators
  Jpos = 1. + braket_pos(J - 1.);
  Jneg = 1. + braket_neg(J - 1.);
  xpos = braket_pos(trbe - 3.);
  xneg = braket_neg(trbe - 3.);

  //compute energies
  Wvolpos = computeWvol(Jpos); Wvolneg = computeWvol(Jneg);
  Wdevpos = computeWinv(xpos); Wdevneg = computeWinv(xneg);

  //compute degradation
  _D[_qp] = (1 - _kdamage[_qp]) * std::pow(1 - _c[_qp], 2. * (1. + _ep[_qp])) + _kdamage[_qp];

  //compose penalized strain energy
  _W[_qp] = (_D[_qp] * Wvolpos) + Wvolneg + (_D[_qp] * Wdevpos) + Wdevneg; 

  //compute history variable

  if ((Wvolpos + Wdevpos) > _Hist_old[_qp]){
    _Hist[_qp] = Wvolpos + Wdevpos;
  }else{
    _Hist[_qp] = _Hist_old[_qp];
  }

  //compute split kirchhoff stress: volumetric part
  Real dWdJpos, dJposdJ;
  Real dWdJneg, dJnegdJ;

  dWdJpos = 0.5 * K * (Jpos - (1. / Jpos)); dJposdJ = (Je > 1) ? 1. : 0.;
  dWdJneg = 0.5 * K * (Jneg - (1. / Jneg)); dJnegdJ = (Je < 1) ? 1. : 0.;

  //define volumetric kirchhoff
  RankTwoTensor tauvol;
  tauvol = _F[_qp].det() * (_D[_qp] * dWdJpos * dJposdJ + dWdJneg * dJnegdJ) * I2;

  //define deviatoric tau and degrade
  RankTwoTensor taudev;
  taudev = _D[_qp] * mu * _be[_qp].deviatoric();

  RankTwoTensor tau = tauvol + taudev;

  _pk1_stress[_qp] = tau * _F[_qp].inverse().transpose();
  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour

  Real ss = std::sqrt(K / _rho[_qp]);
	
  //Compute artificial viscosity term
  Real P_av;
  Real Je_dot;
  Je_dot = ((_Fe[_qp].det() - _Fe_old[_qp].det()) / _dt);
  Real Jtot = _F[_qp].det();

  P_av = _C0 * _rho[_qp] * (Je_dot * std::abs(Je_dot) / std::pow(Je, 2.0)) * std::pow(_Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (Je_dot / Je) * _Le;
  _pk1_stress[_qp] += P_av * I;

  //compute here PK2 stress
  _S[_qp] = _F[_qp].inverse() * _pk1_stress[_qp];

  //cauchy
  _sigma[_qp] = (1. / _Je[_qp]) * tau;
  _Je[_qp] = _Fe[_qp].det();
  _Jp[_qp] = _Fp[_qp].det();
  _Fres[_qp] = _F[_qp] - _Fe[_qp] * _Fp[_qp];
  _E_dot[_qp] = 0.5 * (F_dot.transpose() * F + F.transpose() * F_dot);
  _E[_qp] = 0.5 * (_F[_qp].transpose() * _F[_qp] - I2);

  // Compute the consistent tangent, i.e. the derivative of the PK1 stress w.r.t. the deformation
  // gradient.
  if (_fe_problem.currentlyComputingJacobian())
  {
    RankFourTensor d_tau_d_F = K * detJ * detJ * I.times<i, j, k, l>(Fit) +
                               G * (_d_be_d_F - I.times<i, j, k, l>(I) * _d_be_d_F / 3);
    _pk1_jacobian[_qp] = Fit.times<m, j, i, m, k, l>(d_tau_d_F) - Fit.times<k, j, i, l>(tau * Fit);
  }
  _HS_elastic[_qp] = P_av * _Ee_dot[_qp].trace();
  _HS_plastic[_qp] = _S[_qp].doubleContraction(_Ep_dot[_qp]);
}

Real
NONADLIPIT::computeReferenceResidual(const Real & effective_trial_stress,
                                                              const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  return effective_trial_stress - G * scalar * _be[_qp].trace();
}

Real
NONADLIPIT::computeResidual(const Real & effective_trial_stress,
                                                     const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return effective_trial_stress - G * scalar * _be[_qp].trace() - _H[_qp];
}

Real
NONADLIPIT::computeDerivative(const Real & /*effective_trial_stress*/,
                                                       const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return -G * _be[_qp].trace() - _dH[_qp];
}

void
NONADLIPIT::preStep(const Real & scalar, const Real & R, const Real & J)
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

//define a function to get strain energy
Real
NONADLIPIT::computeWvol(const Real & J)
{
  const Real K = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  //evaluate expression
  //Real Winv = lambda * ((std::pow(I3, 2.) - 1.) / 4.) - ((lambda / 2.) - mu) * std::log(I3);
  Real Wvol = (K / 4.) * (J * J - 1. - 2. * std::log(J));
  return Wvol;
}

Real
NONADLIPIT::computeWinv(const Real & x)
{
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  Real Winv = 0.5 * mu * x;
  return Winv;
}