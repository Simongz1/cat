#include "ADPBXStress.h"

registerMooseObject("mlApp", ADPBXStress);

InputParameters
ADPBXStress::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADComputeStressBase>::validParams();
  params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("Finite Plasticity");
  params.addParam<MaterialPropertyName>("elasticity_tensor", "elasticity_tensor", "The name of the elasticity tensor.");
  params.addRequiredParam<MaterialName>("flow_stress_material","The material defining the flow stress");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredParam<Real>("element_size", "element_size");
  params.addCoupledVar("Yinitial", "Yinitial");
  params.addParam<bool>("euler_angles", true, "euler_angles");
  ///////////////add rule of mixture variables
  params.addRequiredParam<bool>("use_mixture", "use_mixture");
  params.addCoupledVar("fraction_csv", "fraction_csv");
  params.addRequiredParam<Real>("binder_yield", "binder_yield");

  //test with given values for now
  params.addRequiredParam<Real>("binder_bulk", "binder_bulk");
  params.addRequiredParam<Real>("binder_shear", "binder_shear");
  params.addRequiredParam<Real>("poisson_binder", "poisson_binder");
  params.addRequiredCoupledVar("displacements", "displacements");


  //parameters for equations of state
  params.addRequiredParam<Real>("A_unreacted", "A_unreacted");
  params.addRequiredParam<Real>("B_unreacted", "B_unreacted");
  params.addRequiredParam<Real>("R1_unreacted", "R1_unreacted");
  params.addRequiredParam<Real>("R2_unreacted", "R2_unreacted");
  params.addRequiredParam<Real>("omega_unreacted", "omega_unreacted");

  params.addRequiredParam<Real>("A_reacted", "A_reacted");
  params.addRequiredParam<Real>("B_reacted", "B_reacted");
  params.addRequiredParam<Real>("R1_reacted", "R1_reacted");
  params.addRequiredParam<Real>("R2_reacted", "R2_reacted");
  params.addRequiredParam<Real>("omega_reacted", "omega_reacted");
  params.addRequiredCoupledVar("temperature", "temperature");
  return params;
}

ADPBXStress::ADPBXStress(
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

    _be(declareADProperty<RankTwoTensor>("be")),
    _be_old(getMaterialPropertyOldByName<RankTwoTensor>("be")),
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
    _flow_stress_name("yield_mixture"),
    _H(getADMaterialPropertyByName<Real>(_flow_stress_name)),
    _dH(getMaterialPropertyDerivativeByName<Real>(_flow_stress_name, _ep_name)),
    _d2H(getMaterialPropertyDerivativeByName<Real>(_flow_stress_name, _ep_name, _ep_name)),
    /////

    _rho(getADMaterialProperty<Real>("density")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("element_size")),
    _p_unreacted(declareADProperty<Real>("p_unreacted")),
    _p_reacted(declareADProperty<Real>("p_reacted")),
    _p_av(declareADProperty<Real>("p_av")),
    _p_mix(declareADProperty<Real>("p_mix")),

    //for rule of mixture
    _use_mixture(getParam<bool>("use_mixture")),
    _fraction_csv(coupledValue("fraction_csv")),
    _binder_yield(getParam<Real>("binder_yield")),
    _binder_bulk(getParam<Real>("binder_bulk")),
    _binder_shear(getParam<Real>("binder_shear")),

    //other properties
    _inv_Cp(declareADProperty<RankTwoTensor>("inv_Cp")),
    _Cp(declareADProperty<RankTwoTensor>("Cp")),
    _inv_Cp_old(getMaterialPropertyOld<RankTwoTensor>("inv_Cp")),
    _Cp_old(getMaterialPropertyOld<RankTwoTensor>("Cp")),
    _Ce(declareADProperty<RankTwoTensor>("Ce")),
    _Ce_old(getMaterialPropertyOld<RankTwoTensor>("Ce")),
    _Ep_old(getMaterialPropertyOld<RankTwoTensor>("Ep")),
    _Ee_old(getMaterialPropertyOld<RankTwoTensor>("Ee")),
    _S(declareADProperty<RankTwoTensor>("S")),
    _PK1(declareADProperty<RankTwoTensor>("PK1")),
    _HS_elastic(declareADProperty<Real>("HS_elastic")),
    _HS_plastic(declareADProperty<Real>("HS_plastic")),

    //jacobians
    _Je(declareADProperty<Real>("Je")),
    _Jp(declareADProperty<Real>("Jp")),
    _J(declareADProperty<Real>("J")),
    _Je_dot(declareADProperty<Real>("Je_dot")),

    //green lagrange strains
    _E(declareADProperty<RankTwoTensor>("E")),
    _E_dot(declareADProperty<RankTwoTensor>("E_dot")),
    _ndisp(coupledComponents("displacements")),

    //increments
    _strain_increment(getADMaterialProperty<RankTwoTensor>("strain_increment")),
    _rotation_increment(getADMaterialProperty<RankTwoTensor>("rotation_increment")),
    _nu(getParam<Real>("poisson_binder")),
    _ss(declareADProperty<Real>("ss")),
    _Yinitial(adCoupledValue("Yinitial")),

    ////equation of state parameters
    _A_unreacted(getParam<Real>("A_unreacted")),
    _B_unreacted(getParam<Real>("B_unreacted")),
    _R1_unreacted(getParam<Real>("R1_unreacted")),
    _R2_unreacted(getParam<Real>("R2_unreacted")),
    _omega_unreacted(getParam<Real>("omega_unreacted")),

    _A_reacted(getParam<Real>("A_reacted")),
    _B_reacted(getParam<Real>("B_reacted")),
    _R1_reacted(getParam<Real>("R1_reacted")),
    _R2_reacted(getParam<Real>("R2_reacted")),
    _omega_reacted(getParam<Real>("omega_reacted")),
    _temperature(adCoupledValue("temperature")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _pk1_binder(declareADProperty<RankTwoTensor>("pk1_binder"))
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
ADPBXStress::initialSetup()
{
  _flow_stress_material = &getMaterial("flow_stress_material");
}

void
ADPBXStress::initQpStatefulProperties()
{
  ADComputeStressBase::initQpStatefulProperties();
  _be[_qp].setToIdentity();
  _ep[_qp] = 0;
  _Fp[_qp].setToIdentity();
  _Fe[_qp].setToIdentity();
  
  _Cp[_qp].setToIdentity();
  _F[_qp].setToIdentity();
  _Fhat[_qp].setToIdentity();
}

void
ADPBXStress::computeQpStress()
{
  ADRankTwoTensor I2;
  I2.setToIdentity();
  const ADReal G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  const ADReal K = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  const auto I = RankTwoTensor::Identity();

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

  F_incremental = _rotation_increment[_qp] * (_strain_increment[_qp] + I);

  //after deciding method, compute other stuff
  _F[_qp] = F_incremental * _F_old[_qp];
  _Fhat[_qp] = F_incremental;
  Fit = _F[_qp].inverse().transpose();

  ////////////////////////////////////

  // Update configuration
  ADRankTwoTensor f = _Fhat[_qp];
  ADReal Jhat = f.det();
  ADRankTwoTensor f_bar = f / MetaPhysicL::cbrt(Jhat);
  ADReal J = _F[_qp].det();

  // Elastic predictor
  _be[_qp] = f_bar * _be_old[_qp] * f_bar.transpose();

  //here we explicitly use hookes law
  //if rule of mixture
  _mixture_shear = _use_mixture ? 
                            ((_fraction_csv[_qp] * G) + ((1. - _fraction_csv[_qp]) * _binder_shear)) : 
                            G;
  _mixture_bulk = _use_mixture ? 
                            ((_fraction_csv[_qp] * K) + ((1. - _fraction_csv[_qp]) * _binder_bulk)) : 
                            K;
  
  //
  ADRankTwoTensor s = _Yinitial[_qp] * _mixture_shear * _be[_qp].deviatoric();
  ADReal snorm = MetaPhysicL::sqrt(s.doubleContraction(s));
  _Np[_qp] = MooseUtils::absoluteFuzzyEqual(snorm, ADReal(0)) ? std::sqrt(1. / 2.) * I2
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

  //compute plastic strain rate from radial return increment
  _ep_dot[_qp] = delta_ep / _dt;

  ADRankTwoTensor F = _F[_qp];
  ADRankTwoTensor F_bar = MetaPhysicL::pow(J, - 1. / 3.) * F;
  
  _inv_Cp[_qp] = F_bar.inverse() * _be[_qp] * F_bar.inverse().transpose();
  _Cp[_qp] = _inv_Cp[_qp].inverse();

  //get plastic deformation gradient 
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
  _Fe[_qp] = _F[_qp] * _Fp[_qp].inverse();
  _Ce[_qp] = _Fe[_qp].transpose() * _Fe[_qp];
  _Ee[_qp] = 0.5 * (_Ce[_qp] - I2);
  _Ep[_qp] = 0.5 * (_Cp[_qp] - I2);

  ADRankTwoTensor F_dot, Fe_dot, Fp_dot;
  Fp_dot = (1. / _dt) * (_Fp[_qp] - _Fp_old[_qp]);
  Fe_dot = (1. / _dt) * (_Fe[_qp] - _Fe_old[_qp]);
  F_dot = (1. / _dt) * (_F[_qp] - _F_old[_qp]);

  _Ep_dot[_qp] = 0.5 * (Fp_dot.transpose() * _Fp[_qp] + _Fp[_qp].transpose() * Fp_dot);
  _Ee_dot[_qp] = 0.5 * (Fe_dot.transpose() * _Fe[_qp] + _Fe[_qp].transpose() * Fe_dot);

  _Je[_qp] = _Fe[_qp].det();
  _Jp[_qp] = _Fp[_qp].det();
  _J[_qp] = J;

  //compute pressure for unreacted and reacted material
  _p_unreacted[_qp] = - computeJWLPressure(_A_unreacted, _B_unreacted, _R1_unreacted, _R2_unreacted, _omega_unreacted);
  _p_reacted[_qp] = - computeJWLPressure(_A_reacted, _B_reacted, _R1_reacted, _R2_reacted, _omega_reacted);
  _p_mix[_qp] = _Yinitial[_qp] * _p_unreacted[_qp] + (1. - _Yinitial[_qp]) * _p_reacted[_qp];

  ///BINDER MODEL///
  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const Real beta = (2 * _nu) / (1 - 2 * _nu);

  //binder stress calculation
  ADRankTwoTensor FJF = _Fe[_qp] - MetaPhysicL::pow(_Fe[_qp].det(), - beta) * _Fe[_qp].inverse().transpose();
  ADRankTwoTensor pk1_binder = _binder_shear * FJF * _Fp[_qp].inverse().transpose();
  ADRankTwoTensor tau_binder = pk1_binder * _F[_qp].transpose();

  //output binder stress
  _pk1_binder[_qp] = pk1_binder;

  //RDX tau
  ADRankTwoTensor tau_RDX = _p_mix[_qp] * I + s;
                          
  ADRankTwoTensor tau_total;
  
  if (_use_mixture){
    tau_total = _fraction_csv[_qp] * tau_RDX + (1. - _fraction_csv[_qp]) * tau_binder;
  }else{
    tau_total = tau_RDX;
  }

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour
  ADReal ss = MetaPhysicL::sqrt(_mixture_bulk / _rho[_qp]);
  _ss[_qp] = ss;

  //compute artificial viscosity
  _p_av[_qp] = computeAVPressure();
	
  //cauchy
  _stress[_qp] = (1. / _J[_qp]) * tau_total + _p_av[_qp] * I;
  _PK1[_qp] = _J[_qp] * _stress[_qp] * _F[_qp].inverse().transpose();

  //PK1 total
  //ADRankTwoTensor pk1_total = _Je[_qp] * _stress[_qp] * _Fe[_qp].inverse().transpose();

  //compute here PK2 stress
  _S[_qp] = _F[_qp].inverse() * _PK1[_qp];

  
  _E_dot[_qp] = 0.5 * (F_dot.transpose() * F + F.transpose() * F_dot);
  _E[_qp] = 0.5 * (_F[_qp].transpose() * _F[_qp] - I2);

  _HS_elastic[_qp] = _p_av[_qp] * _Ee_dot[_qp].trace();
  _HS_plastic[_qp] = _S[_qp].doubleContraction(_Ep_dot[_qp]);
}

Real
ADPBXStress::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                              const ADReal & scalar)
{
  const ADReal G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  return MetaPhysicL::raw_value(effective_trial_stress - _mixture_shear * scalar * _be[_qp].trace());
}

ADReal
ADPBXStress::computeResidual(const ADReal & effective_trial_stress,
                                                     const ADReal & scalar)
{
  const ADReal G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return (effective_trial_stress - _mixture_shear * scalar * _be[_qp].trace() - _H[_qp]);
}

ADReal
ADPBXStress::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                       const ADReal & scalar)
{
  const ADReal G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return (- _mixture_shear * _be[_qp].trace() - _dH[_qp]);
}

ADReal
ADPBXStress::computeJWLPressure(const Real & A, 
                                      const Real & B,
                                      const Real & R1, 
                                      const Real & R2, 
                                      const Real & omega)
{
  ADReal p = A * MetaPhysicL::exp( - R1 * _Je[_qp] ) + B * MetaPhysicL::exp( - R2 * _Je[_qp]);
  p += omega * _rho[_qp] * _cv[_qp] * _temperature[_qp] / _Je[_qp];
  return p;
}

ADReal
ADPBXStress::computeAVPressure()
{
  ADReal P_av;
  ADReal J = _F[_qp].det();
  ADReal J_dot = ((_F[_qp].det() - _F_old[_qp].det()) / _dt);

  ADReal Le = ADReal(_Le);

  P_av = _C0 * _rho[_qp] * (J_dot * MetaPhysicL::abs(J_dot) / MetaPhysicL::pow(J, 2.0)) * MetaPhysicL::pow(Le, 2.0);
  P_av += _C1 * _rho[_qp] * _ss[_qp] * (J_dot / J) * Le;
  return P_av;
}
