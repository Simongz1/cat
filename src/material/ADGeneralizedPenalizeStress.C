#include "ADGeneralizedPenalizeStress.h"

registerMooseObject("mistApp", ADGeneralizedPenalizeStress);

InputParameters
ADGeneralizedPenalizeStress::validParams()
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

  //test
  params.addParam<MaterialPropertyName>("D_name", "D_name");
  params.addParam<VariableName>("c_name", "c_name");

  params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "kappa property name");
  params.addParam<MaterialPropertyName>("mobility_name", "L", "mobility property name");
  params.addParam<MaterialPropertyName>("elastic_energy_name", "elastic_energy", "elastic energy property name");
  params.addRequiredParam<Real>("ep_ref", "refference effective plastic strain");
  params.addRequiredCoupledVar("h_min", "h_min");

  //EOS and chemistry stuff
  params.addCoupledVar("Y_unreacted", "Y_unreacted");
  return params;
}

ADGeneralizedPenalizeStress::ADGeneralizedPenalizeStress(
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
    _c(adCoupledValue("c")),
    _l(getParam<Real>("l")),
    _kappa(declareProperty<Real>(getParam<MaterialPropertyName>("kappa_name"))),
    _L(declareProperty<Real>(getParam<MaterialPropertyName>("mobility_name"))),
    _gc(getParam<Real>("gc")),
    _gcprop(coupledValue("gcprop")),
    _visco(getParam<Real>("visco")),
    _kdamage(getParam<Real>("kdamage")),

    //penalization function sent to the AD system
    _D(declareADProperty<Real>("D")),

    _Hist(declareProperty<Real>("Hist")),
    _Hist_old(getMaterialPropertyOld<Real>("Hist")),
    _elastic_energy(declareADProperty<Real>("elastic_energy")),
    _delastic_energydc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("elastic_energy_name"), coupledName("c", 0))),
    _d2elastic_energyd2c(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("elastic_energy_name"), coupledName("c", 0), coupledName("c", 0))),
    _dstress_dc(declarePropertyDerivative<Real>(_base_name + "stress", coupledName("c", 0))),
    _pk1_pos(declareProperty<RankTwoTensor>("pk1_pos")),
    _pk1_neg(declareProperty<RankTwoTensor>("pk1_neg")),
    _W(declareADProperty<Real>("W")),
    _Wpos(declareADProperty<Real>("Wpos")),
    _Wneg(declareADProperty<Real>("Wneg")),
    //invariants for debugging
    _elastic_energy_total(declareADProperty<Real>("elastic_energy_total")),
    _ep_ref(getParam<Real>("ep_ref")),
    _Cp_bar(declareProperty<RankTwoTensor>("Cp_bar")),
    _Cp(declareProperty<RankTwoTensor>("Cp")),
    _Cp_bar_old(getMaterialPropertyOld<RankTwoTensor>("Cp_bar")),
    _Cp_old(getMaterialPropertyOld<RankTwoTensor>("Cp")),
    _Ep_old(getMaterialPropertyOld<RankTwoTensor>("Ep")),
    _Ee_old(getMaterialPropertyOld<RankTwoTensor>("Ee")),
    _F_computed(declareProperty<RankTwoTensor>("F_computed")),
    _S(declareProperty<RankTwoTensor>("S")),
    _HS_elastic(declareADProperty<Real>("HS_elastic")),
    _HS_plastic(declareADProperty<Real>("HS_plastic")),
    _C_computed(declareProperty<RankTwoTensor>("C_computed")),

    //EOS and chemistry stuff
    _pressure_total(getADMaterialProperty<Real>("pressure_total")),
    _Y_unreacted(adCoupledValue("Y_unreacted"))
{
}

void
ADGeneralizedPenalizeStress::initialSetup()
{
  _flow_stress_material = &getMaterial("flow_stress_material");

  // Enforce isotropic elastic tensor
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("ADGeneralizedPenalizeStress requires an isotropic elasticity tensor");
}

void
ADGeneralizedPenalizeStress::initQpStatefulProperties()
{
  ComputeLagrangianStressPK1::initQpStatefulProperties();
  _be[_qp].setToIdentity();
  _ep[_qp] = 0;
  _Fp[_qp].setToIdentity(); //stateful plastic deformation gradient
  _Fe[_qp].setToIdentity();
  _Cp_bar[_qp].zero();
  _Cp[_qp].zero();
}

void
ADGeneralizedPenalizeStress::computeQpPK1Stress()
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

  //obtain inverse plastic volume preserving C tensor

  //compute F_bar at n+1
  //RankTwoTensor F = f * _F_old[_qp].inverse(); //equivalent to computing F[n+1] = f[n+1]F[n].inverse();
  RankTwoTensor F = _deformation_gradient[_qp];
  _F_computed[_qp] = F;
  RankTwoTensor F_bar = std::pow(F.det(), - 1. / 3.) * F;

  //use identity to get volume preserving C^p^-1

  _Cp_bar[_qp] = F_bar.inverse() * _be[_qp] * F_bar.inverse().transpose();
  _Cp_bar[_qp] = _Cp_bar[_qp].inverse();
  _Cp[_qp] = _Cp_bar[_qp];

  //use this to compute plastic strain

  _Ep[_qp] = 0.5 * (_Cp[_qp] - I);
  _Ep_dot[_qp] = (1. / _dt) * (_Ep[_qp] - _Ep_old[_qp]); //backward scheme
  
  RankTwoTensor be_total = F * _Cp[_qp].inverse() * F.transpose();
  _Ee[_qp] = 0.5 * (be_total.transpose() - I);
  _Ee_dot[_qp] = (1. / _dt) * (_Ee[_qp] - _Ee_old[_qp]);

  ///invariants for elastic energy calculation

  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);

  Real I1;
  Real I3;
  RankTwoTensor B = F * F.transpose();
  RankTwoTensor C = F.transpose() * F;

  _C_computed[_qp] = C;

  //elastic energy stuff for fracture
  _kappa[_qp] = _gcprop[_qp] * _l;
  _L[_qp] = 1. / (_gcprop[_qp] * _visco);

  //assign degradation derivatives
  ADReal S = 1. - _c[_qp]; //degradation value
  ADReal e_norm = _ep[_qp] / _ep_ref; //normalized plastic strain
  ADReal D = std::pow(S, 2 * e_norm) + _kdamage;
  //store
  _D[_qp] = D;
  //
  ADReal dDdc = - 2. * e_norm * std::pow(S, (2.* e_norm) - 1.); //derivative of degradation w.r.t. c
  ADReal d2Dd2c = 2. * e_norm * (2. * e_norm - 1.) * std::pow(S, (2. * e_norm) - 2.);

  Real elastic_energy = lambda * ((std::pow(F.det(), 2.) - 1.) / 4.) - ((lambda / 2.) - mu) * std::log(F.det()) + 0.5 * mu * (C.trace() - 3.);
  _elastic_energy_total[_qp] = elastic_energy; //undifferentiated elastic energy directly from cauchy tensor

  //test: penalize shear 
  _pk1_stress[_qp] = _pressure_total[_qp].value() * I + _D[_qp].value() * _Y_unreacted[_qp].value() * G * _be[_qp].deviatoric();
  _S[_qp] = _F_computed[_qp].inverse() * _pk1_stress[_qp];

  /////////////COMPUTE HEAT SOURCES HERE

  _HS_plastic[_qp] = 0.5 * _Y_unreacted[_qp] *  _S[_qp].doubleContraction(_Ep_dot[_qp]);

  //test: decompose strain energy based in sign

  _Wpos[_qp] = std::max(_elastic_energy_total[_qp], 0.);
  _Wneg[_qp] = _elastic_energy_total[_qp] - _Wpos[_qp];

  //penalize positive
  _W[_qp] = _D[_qp] * _Wpos[_qp] + _Wneg[_qp];

  if (_Wpos[_qp].value() > _Hist_old[_qp]){
    _Hist[_qp] = _Wpos[_qp].value();
  }else{
    _Hist[_qp] = _Hist_old[_qp];
  }

  //damage stuff

  _dstress_dc[_qp] = D.value();
  _elastic_energy[_qp] = (D * _Hist[_qp]) + (_gcprop[_qp] * std::pow(_c[_qp], 2.) / (2. * _l)); //sum of penalized elastic energy plus fractured new surface energy
  _delastic_energydc[_qp] = (dDdc.value() * _Hist[_qp]) + (_gcprop[_qp] * _c[_qp].value() / _l);
  _d2elastic_energyd2c[_qp] = (d2Dd2c.value() * _Hist[_qp]) + (_gcprop[_qp] / _l);
  
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
  _pk1_stress[_qp] += P_av * I;

  //define tau

  RankTwoTensor tau = _pk1_stress[_qp] * _F_computed[_qp].transpose();
  // Compute the consistent tangent, i.e. the derivative of the PK1 stress w.r.t. the deformation
  // gradient.
  if (_fe_problem.currentlyComputingJacobian())
  {
    RankFourTensor d_tau_d_F = K * detJ * detJ * I.times<i, j, k, l>(Fit) +
                               G * (_d_be_d_F - I.times<i, j, k, l>(I) * _d_be_d_F / 3);
    _pk1_jacobian[_qp] = Fit.times<m, j, i, m, k, l>(d_tau_d_F) - Fit.times<k, j, i, l>(tau * Fit);
  }
  _HS_elastic[_qp] = 0.5 * P_av * _Ee_dot[_qp].trace();
}

Real
ADGeneralizedPenalizeStress::computeReferenceResidual(const Real & effective_trial_stress,
                                                              const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  return effective_trial_stress - G * scalar * _be[_qp].trace();
}

Real
ADGeneralizedPenalizeStress::computeResidual(const Real & effective_trial_stress,
                                                     const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return effective_trial_stress - G * scalar * _be[_qp].trace() - _H[_qp];
}

Real
ADGeneralizedPenalizeStress::computeDerivative(const Real & /*effective_trial_stress*/,
                                                       const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  return -G * _be[_qp].trace() - _dH[_qp];
}

void
ADGeneralizedPenalizeStress::preStep(const Real & scalar, const Real & R, const Real & J)
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