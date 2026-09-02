#include "ADPBXStress.h"

registerMooseObject("mlApp", ADPBXStress);

InputParameters
ADPBXStress::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADComputeStressBase>::validParams();
  params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("Object that computes finite strain elasto-plastic response.");
  params.addRequiredParam<MaterialName>("flow_stress_material","The material defining the flow stress");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredParam<Real>("element_size", "element_size");
  params.addCoupledVar("Yinitial", "Yinitial");
  params.addParam<bool>("euler_angles", true, "euler_angles");
  ///////////////add rule of mixture variables
  params.addRequiredParam<bool>("use_mixture", "use_mixture");
  params.addCoupledVar("fraction_csv", "fraction_csv");

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
    _F(declareADProperty<RankTwoTensor>("deformation_gradient")),
    _F_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _Fhat(declareADProperty<RankTwoTensor>("Fhat")),

    //obtain bulk and shear for RDX and binder
    _binder_bulk(getADMaterialProperty<Real>("binder_bulk")),
    _binder_poisson(getADMaterialProperty<Real>("binder_poisson")),
    _binder_yield(getADMaterialProperty<Real>("binder_yield")),
    _sigma_binder(declareADProperty<RankTwoTensor>("sigma_binder")),

    _RDX_bulk(getADMaterialProperty<Real>("RDX_bulk")),
    _RDX_poisson(getADMaterialProperty<Real>("RDX_poisson")),
    _sigma_RDX(declareADProperty<RankTwoTensor>("sigma_RDX")),

    //properties for plasticity
    _ep_name("ep"),
    _ep(declareADProperty<Real>(_ep_name)),
    _ep_old(getMaterialPropertyOldByName<Real>(_ep_name)),
    _ep_dot(declareADProperty<Real>("ep_dot")),

    _be_bar(declareADProperty<RankTwoTensor>("be_bar")),
    _be_bar_old(getMaterialPropertyOldByName<RankTwoTensor>("be_bar")),
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

    //other properties
    _Cp(declareADProperty<RankTwoTensor>("Cp")),

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

    //green lagrange strains
    _E(declareADProperty<RankTwoTensor>("E")),
    _E_dot(declareADProperty<RankTwoTensor>("E_dot")),

    //increments
    _strain_increment(getADMaterialProperty<RankTwoTensor>("strain_increment")),
    _rotation_increment(getADMaterialProperty<RankTwoTensor>("rotation_increment")),

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
    _cv(getADMaterialProperty<Real>("specific_heat"))
{}

void
ADPBXStress::initialSetup()
{
  _flow_stress_material = &getMaterial("flow_stress_material");
}

void
ADPBXStress::initQpStatefulProperties()
{
  ADComputeStressBase::initQpStatefulProperties();
  _be_bar[_qp].setToIdentity();
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
  //compute the shear modulus from bulk and poisson
  const ADReal RDX_shear = (3.0 * _RDX_bulk[_qp]) * (1.0 - 2.0 * _RDX_poisson[_qp]) / (2.0 + 2.0 * _RDX_poisson[_qp]);
  const ADReal binder_shear = (3.0 * _binder_bulk[_qp]) * (1.0 - 2.0 * _binder_poisson[_qp]) / (2.0 + 2.0 * _binder_poisson[_qp]);

  ADRankTwoTensor I2;
  I2.setToIdentity();
  const auto I = RankTwoTensor::Identity();
  //compute AD version of the incremental deformation gradient
  //form a tensor with rows occupied by displacements

  ADRankTwoTensor F_incremental;
  F_incremental = _rotation_increment[_qp] * (_strain_increment[_qp] + I);

  //after deciding method, compute other stuff
  _F[_qp] = F_incremental * _F_old[_qp];
  _Fhat[_qp] = F_incremental;
  ////////////////////////////////////

  //update configuration
  ADRankTwoTensor f = _Fhat[_qp];
  ADReal Jhat = f.det();
  ADRankTwoTensor f_bar = f / MetaPhysicL::cbrt(Jhat);
  ADReal J = _F[_qp].det();

  //elastic predictor
  _be_bar[_qp] = f_bar * _be_bar_old[_qp] * f_bar.transpose();

  //compute INITIAL mixture properties as function of static fractions
  _mixture_shear = _use_mixture ? 
                            _fraction_csv[_qp] * RDX_shear + (1 - _fraction_csv[_qp]) * binder_shear : 
                            RDX_shear;
  _mixture_bulk = _use_mixture ? 
                            _fraction_csv[_qp] * _RDX_bulk[_qp] + (1 - _fraction_csv[_qp]) * _binder_bulk[_qp] : 
                            _RDX_bulk[_qp];
  
  ADRankTwoTensor s = _Yinitial[_qp] * _mixture_shear * _be_bar[_qp].deviatoric();
  ADReal snorm = MetaPhysicL::sqrt(s.doubleContraction(s));
  _Np[_qp] = MooseUtils::absoluteFuzzyEqual(snorm, ADReal(0)) ? std::sqrt(1. / 2.) * I2
                                                         : std::sqrt(3. / 2.) * s / snorm;
  ADReal s_eff = s.doubleContraction(_Np[_qp]);

  // Check for plastic loading and do return mapping
  ADReal delta_ep = 0;

  //run radial return 
  if (MetaPhysicL::raw_value(computeResidual(s_eff, 0)) > 0)
  {
    returnMappingSolve(s_eff, delta_ep, _console);
  }

  // Update intermediate and current configurations
  _ep[_qp] = _ep_old[_qp] + delta_ep;
  _be_bar[_qp] -= 2. / 3. * delta_ep * _be_bar[_qp].trace() * _Np[_qp];

  // Recompute the effective deviatoric stress after the plastic correction.
  // The return mapping is formulated for the homogenized material, so this
  // must use the same mixture shear modulus as the trial stress and residual.
  ADRankTwoTensor s_corrected =
      _Yinitial[_qp] * _mixture_shear * _be_bar[_qp].deviatoric();
    
  //compute plastic strain rate from radial return increment
  _ep_dot[_qp] = delta_ep / _dt;

  //here, we recover the actual be, not the volume preserving part
  ADRankTwoTensor be = MetaPhysicL::pow(J, 2.0 / 3.0) * _be_bar[_qp];
  _Cp[_qp] = _F[_qp].transpose() * be.inverse() * _F[_qp];
  
  //ASSUMING that there is no plastic rotation
  //we approximate the plastic deformation gradient as the sqrt of Cp
  //here we obtain the Ce tensor using a polar decomposition
  //get the symmetric be
  ADRankTwoTensor Cp_sym = 0.5 * (_Cp[_qp] + _Cp[_qp].transpose());

  //get plastic deformation gradient 
  ADRankTwoTensor Q;
  std::vector<ADReal> lam(3);

  //obtain eigenvalues and eigenvectors
  Cp_sym.symmetricEigenvaluesEigenvectors(lam, Q);

  //obtain the diagonal tensor with eigenvalues as principal diagonal
  ADRankTwoTensor sqrt_diag; sqrt_diag.zero();

  //populate the square root of the diagonal tensor
  for (unsigned int i = 0; i < 3; ++i){
    sqrt_diag(i,i) = MetaPhysicL::sqrt(std::max(lam[i], ADReal(1e-12)));
  }

  //reconstruct plastic deformation gradient
  _Fp[_qp] = Q * sqrt_diag * Q.transpose();
  //recover elastic deformation gradient
  _Fe[_qp] = _F[_qp] * _Fp[_qp].inverse();

  //compute Lagrangian strains
  _Ee[_qp] = 0.5 * (_Fe[_qp].transpose() * _Fe[_qp] - I2);
  _Ep[_qp] = 0.5 * (_Cp[_qp] - I2);

  ADRankTwoTensor F_dot, Fe_dot, Fp_dot;
  Fp_dot = (1. / _dt) * (_Fp[_qp] - _Fp_old[_qp]);
  Fe_dot = (1. / _dt) * (_Fe[_qp] - _Fe_old[_qp]);
  F_dot = (1. / _dt) * (_F[_qp] - _F_old[_qp]);

  _Ep_dot[_qp] = 0.5 * (Fp_dot.transpose() * _Fp[_qp] + _Fp[_qp].transpose() * Fp_dot);
  _Ee_dot[_qp] = 0.5 * (Fe_dot.transpose() * _Fe[_qp] + _Fe[_qp].transpose() * Fe_dot);

  _Je[_qp] = _Fe[_qp].det();
  _Jp[_qp] = _Fp[_qp].det(); //this is only for checking consistency
  _J[_qp] = J;

  //compute pressure for unreacted and reacted material
  _p_unreacted[_qp] = - computeJWLPressure(_A_unreacted, _B_unreacted, _R1_unreacted, _R2_unreacted, _omega_unreacted);
  _p_reacted[_qp] = - computeJWLPressure(_A_reacted, _B_reacted, _R1_reacted, _R2_reacted, _omega_reacted);
  _p_mix[_qp] = _Yinitial[_qp] * _p_unreacted[_qp] + (1. - _Yinitial[_qp]) * _p_reacted[_qp];

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour
  ADReal ss = MetaPhysicL::sqrt(_mixture_bulk + ((4. / 3.) * _mixture_shear) / (_rho[_qp] / _J[_qp]) );
  _ss[_qp] = ss;

  //compute artificial viscosity
  _p_av[_qp] = computeAVPressure();

  //compute the binder stress using a compressible neo-hookean model
  _sigma_binder[_qp] = _binder_bulk[_qp] * (_J[_qp] - 1.0) * I2;
  _sigma_binder[_qp] += (1.0 / _J[_qp]) * binder_shear * _be_bar[_qp].deviatoric();

  //declare stress tensor for the RDX phase ONLY
  _sigma_RDX[_qp] = s_corrected + _p_mix[_qp] * I;
	
  // The deviatoric response is already homogenized through _mixture_shear.
  // Apply the phase rule of mixtures only to the volumetric response; mixing
  // the complete phase stresses here would weight the shear response twice.
  if (_use_mixture)
  {
    const ADRankTwoTensor sigma_binder_vol =
        _binder_bulk[_qp] * (_J[_qp] - 1.0) * I2;
    const ADRankTwoTensor sigma_RDX_vol = _p_mix[_qp] * I;
    _stress[_qp] = s_corrected +
                   (1.0 - _fraction_csv[_qp]) * sigma_binder_vol +
                   _fraction_csv[_qp] * sigma_RDX_vol;
  }
  else
    _stress[_qp] = _sigma_RDX[_qp];

  _stress[_qp] += _p_av[_qp] * I;
  _PK1[_qp] = _J[_qp] * _stress[_qp] * _F[_qp].inverse().transpose();

  //PK1 total
  //ADRankTwoTensor pk1_total = _Je[_qp] * _stress[_qp] * _Fe[_qp].inverse().transpose();

  //compute here PK2 stress. This is needed for computing work conjugate quantities.
  _S[_qp] = _F[_qp].inverse() * _PK1[_qp];
  
  //compute TOTAL LAGRANGIAN strain and strain rate.
  _E_dot[_qp] = 0.5 * (F_dot.transpose() * _F[_qp] + _F[_qp].transpose() * F_dot);
  _E[_qp] = 0.5 * (_F[_qp].transpose() * _F[_qp] - I2);
  
  //compute derived heat sources from these quantities.
  _HS_elastic[_qp] = _p_av[_qp] * _Ee_dot[_qp].trace();
  _HS_plastic[_qp] = _S[_qp].doubleContraction(_Ep_dot[_qp]);
}

Real
ADPBXStress::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                              const ADReal & scalar)
{
  return MetaPhysicL::raw_value(
      effective_trial_stress - _Yinitial[_qp] * _mixture_shear * scalar * _be_bar[_qp].trace());
}

ADReal
ADPBXStress::computeResidual(const ADReal & effective_trial_stress,
                                                     const ADReal & scalar)
{
  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);
  
  //here we require an exact coupling between Y and the defradation of the yield surface
  //we coupld possibly define degradation functions as in phase field
  return (effective_trial_stress - _Yinitial[_qp] * _mixture_shear * scalar * _be_bar[_qp].trace() - _H[_qp]);
}

ADReal
ADPBXStress::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                       const ADReal & scalar)
{
  //update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return (-_Yinitial[_qp] * _mixture_shear * _be_bar[_qp].trace() - _dH[_qp]);
}

ADReal
ADPBXStress::computeJWLPressure(const Real & A, 
                                      const Real & B,
                                      const Real & R1, 
                                      const Real & R2, 
                                      const Real & omega)
{
  ADReal p = A * MetaPhysicL::exp( - R1 * _Je[_qp] ) + B * MetaPhysicL::exp( - R2 * _Je[_qp]);
  p += omega * _rho[_qp] * _cv[_qp] * (_temperature[_qp] - 300.0) / _Je[_qp];
  p -= A * std::exp(- R1) + B * std::exp(- R2);
  return p;
}

ADReal
ADPBXStress::computeAVPressure()
{ 
  //immediately check for compression
  ADReal J_dot = ((_F[_qp].det() - _F_old[_qp].det()) / _dt);

  if (MetaPhysicL::raw_value(J_dot) >= 0.0){
    //this means expansion, then immediately return 0
    return 0.0;
  }
  
  //else, compute the normal expression 
  else{
    ADReal P_av;
    ADReal J = _F[_qp].det();
    ADReal Le = ADReal(_Le);

    P_av = _C0 * _rho[_qp] * (J_dot * MetaPhysicL::abs(J_dot) / MetaPhysicL::pow(J, 2.0)) * MetaPhysicL::pow(Le, 2.0);
    P_av += _C1 * _rho[_qp] * _ss[_qp] * (J_dot / J) * Le;
    return P_av;
  }
}
//END