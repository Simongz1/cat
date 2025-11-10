#include "ADAnisoStress.h"

registerMooseObject("mistApp", ADAnisoStress);

//test: currently changing Yfinal to retrieve ADVariableValue

InputParameters
ADAnisoStress::validParams()
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
  params.addCoupledVar("Yinitial", "Yinitial");
  params.addParam<bool>("euler_angles", true, "euler_angles");
  params.addCoupledVar("euler1", "euler1");
  params.addCoupledVar("euler2", "euler2");
  params.addCoupledVar("euler3", "euler3");
  return params;
}

ADAnisoStress::ADAnisoStress(
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
    //_pressure_av(declareProperty<RankTwoTensor>("pressure_av")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _pressure_total(getADMaterialProperty<Real>("pressure_total")),

    _plastic_strain(declareProperty<RankTwoTensor>("plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("plastic_strain")),
    _plastic_strain_rate(declareProperty<RankTwoTensor>("plastic_strain_rate")),
    _pnorm(declareProperty<Real>("pnorm")),
    _hsp(declareProperty<Real>("hsp")),
    _cauchy_stress(getMaterialProperty<RankTwoTensor>("cauchy_stress")),
    _Yinitial(adCoupledValue("Yinitial")),

    //get euler fractions
    _euler1(coupledValue("euler1")),
    _euler2(coupledValue("euler2")),
    _euler3(coupledValue("euler3"))
    
{
}

void
ADAnisoStress::initialSetup()
{
  _flow_stress_material = &getMaterial("flow_stress_material");

  // Enforce isotropic elastic tensor
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseWarning("ADAnisoStress requires an isotropic elasticity tensor");
}

void
ADAnisoStress::initQpStatefulProperties()
{
  ComputeLagrangianStressPK1::initQpStatefulProperties();
  _be[_qp].setToIdentity();
  _ep[_qp] = 0;
  _plastic_strain[_qp].setToIdentity();
  _Fp[_qp].setToIdentity(); //stateful plastic deformation gradient
  _Fe[_qp].setToIdentity();
}

void
ADAnisoStress::computeQpPK1Stress()
{
  //compute the rotated elasticity tensor
  const RankTwoTensor R = getRotationMatrix(_euler1[_qp] * 2. * M_PI, _euler2[_qp] * 2. * M_PI, _euler3[_qp] * 2. * M_PI);
  _rotated_elasticity_tensor = _elasticity_tensor[_qp];
  _rotated_elasticity_tensor.rotate(R.transpose());

  usingTensorIndices(i, j, k, l, m);
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_rotated_elasticity_tensor);
  const Real K = ElasticityTensorTools::getIsotropicBulkModulus(_rotated_elasticity_tensor);
  const auto I = RankTwoTensor::Identity();
  const auto Fit = _F[_qp].inverse().transpose();
  const auto detJ = _F[_qp].det();

  // Update configuration
  RankTwoTensor f = _inv_df[_qp].inverse();
  RankTwoTensor f_bar = f / std::cbrt(f.det());

  // Elastic predictor
  _be[_qp] = f_bar * _be_old[_qp] * f_bar.transpose();

  //here we explicitly use hookes law
  RankTwoTensor s = _rotated_elasticity_tensor * _be[_qp].deviatoric();
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

  //compute shear and kirchhoff stresses
  s = _Yinitial[_qp].value() * (_rotated_elasticity_tensor * _be[_qp].deviatoric()); //changed to Yinitial. Only unreacted phase can affort shear stress
  RankTwoTensor tau = _pressure_total[_qp].value() * I + s;
  _pk1_stress[_qp] = tau * Fit;

  //compute lagrangian plastic strain from plastic part of the deformation gradient
  _plastic_strain[_qp] = (1. / 2.) * (_Fp[_qp].transpose() * _Fp[_qp] - I);
  _pnorm[_qp] = std::sqrt(2. / 3.) * std::sqrt(_plastic_strain[_qp].doubleContraction(_plastic_strain[_qp]));
  _plastic_strain_rate[_qp] = (_plastic_strain[_qp] - _plastic_strain_old[_qp]) / _dt;
  _hsp[_qp] = _cauchy_stress[_qp].doubleContraction(_plastic_strain_rate[_qp]);

  // Compute the consistent tangent, i.e. the derivative of the PK1 stress w.r.t. the deformation
  // gradient.
  if (_fe_problem.currentlyComputingJacobian())
  {
    RankFourTensor d_tau_d_F = K * detJ * detJ * I.times<i, j, k, l>(Fit) +
                               G * (_d_be_d_F - I.times<i, j, k, l>(I) * _d_be_d_F / 3);
    _pk1_jacobian[_qp] = Fit.times<m, j, i, m, k, l>(d_tau_d_F) - Fit.times<k, j, i, l>(tau * Fit);
  }

  //compute artificial viscosity after stress computations

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
}

Real
ADAnisoStress::computeReferenceResidual(const Real & effective_trial_stress,
                                                              const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_rotated_elasticity_tensor);

  //same

  RankTwoTensor CN = _rotated_elasticity_tensor * _Np[_qp];
  Real NCN = _Np[_qp].doubleContraction(CN);

  return effective_trial_stress - scalar * NCN;
}

Real
ADAnisoStress::computeResidual(const Real & effective_trial_stress,
                                                     const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_rotated_elasticity_tensor);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  //modification for anisotropic plasticity formulation

  RankTwoTensor CN = _rotated_elasticity_tensor * _Np[_qp];
  Real NCN = _Np[_qp].doubleContraction(CN);

  return effective_trial_stress - scalar * NCN - _H[_qp];
}

Real
ADAnisoStress::computeDerivative(const Real & /*effective_trial_stress*/,
                                                       const Real & scalar)
{
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_rotated_elasticity_tensor);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  //same as above

  RankTwoTensor CN = _rotated_elasticity_tensor * _Np[_qp];
  Real NCN = _Np[_qp].doubleContraction(CN);

  return - NCN - _dH[_qp];
}

void
ADAnisoStress::preStep(const Real & scalar, const Real & R, const Real & J)
{
  if (!_fe_problem.currentlyComputingJacobian())
    return;

  const auto I = RankTwoTensor::Identity();
  const Real G = ElasticityTensorTools::getIsotropicShearModulus(_rotated_elasticity_tensor);

  // Update the flow stress
  _ep[_qp] = _ep_old[_qp] + scalar;
  _flow_stress_material->computePropertiesAtQp(_qp);

  _d_R_d_betr =
      G * _Np[_qp] - G * scalar * I - (G * _be[_qp].trace() + _dH[_qp]) * _d_deltaep_d_betr;
  _d_J_d_betr = -G * I - _d2H[_qp] * _d_deltaep_d_betr;
  _d_deltaep_d_betr += -1 / J * _d_R_d_betr + R / J / J * _d_J_d_betr;
}

//write function to generate rotation matrix
RankTwoTensor 
ADAnisoStress::getRotationMatrix(const Real euler1, const Real euler2, const Real euler3){
  //initially, we form the rotation about each angle
  RankTwoTensor R1;
  RankTwoTensor R2;
  RankTwoTensor R3;

  //rotation about x axis
  R1(0,0) = std::cos(euler1);
  R1(0,1) = - std::sin(euler1);
  R1(0,2) = 0.;

  R1(1,0) = - R1(0,1);
  R1(1,1) = R1(0,0);
  R1(1,2) = 0.;

  R1(2,0) = 0.;
  R1(2,1) = 0.;
  R1(2,2) = 1.;

  //rotation about y axis
  R2(0,0) = std::cos(euler2);
  R2(0,1) = 0.;
  R2(0,2) = std::sin(euler2);

  R2(1,0) = 0.;
  R2(1,1) = 1.;
  R2(1,2) = 0.;

  R2(2,0) = - R2(0,2);
  R2(2,1) = 0.;
  R2(2,2) = R2(0,0);
  
  //rotation about z axis
  R3(0,0) = 1.;
  R3(0,1) = 0.;
  R3(0,2) = 0.;

  R3(1,0) = 0.;
  R3(1,1) = std::cos(euler3);
  R3(1,2) = - std::sin(euler3);

  R3(2,0) = 0.;
  R3(2,1) = - R3(1,2);
  R3(2,2) = R3(1,1);

  //now return the product
  return R3 * R1 * R2;
} 