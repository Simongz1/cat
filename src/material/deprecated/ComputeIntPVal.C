#include "ComputeIntPVal.h"

registerMooseObject("TensorMechanicsApp", ComputeIntPVal);

InputParameters
ComputeIntPVal::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Standard compute Mie Gruneisen Pressure with JWL pressure for reacted material. Also computes artificial viscosity contribution");
  params.addRequiredParam<Real>("Gamma", "Grunseisen Parameter");
  params.addRequiredParam<Real>("T_ref", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("element_size", "element_size");
  params.addRequiredParam<Real>("slope", "slope");
  params.addRequiredCoupledVar("Y_final", "final products mass fraction");
  params.addRequiredParam<Real>("A1", "JWL reacted EOS parameter A1");
  params.addRequiredParam<Real>("R1", "JWL reacted EOS parameter B1");
  params.addRequiredParam<Real>("A2", "JWL reacted EOS parameter A2");
  params.addRequiredParam<Real>("R2", "JWL reacted EOS parameter B2");
  params.addRequiredParam<Real>("omega", "JWL reacted EOS parameter omega");
  params.addRequiredParam<Real>("flag_threshold", "pressure flag threshold value");
  return params;
}

ComputeIntPVal::ComputeIntPVal(
    const InputParameters & parameters)
  : Material(parameters),
    _Gamma(getParam<Real>("Gamma")),
    _T_ref(getParam<Real>("T_ref")),
    _rho(getMaterialProperty<Real>("density")),
    _Cv(getMaterialProperty<Real>("specific_heat")),
    _T(coupledValue("temperature")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>("mechanical_strain")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _Le(getParam<Real>("element_size")),
    _pressure_mg(declareProperty<Real>("pressure_mg")),
    _pressure_JWL(declareProperty<Real>("pressure_JWL")),
    _pressure_total(declareProperty<Real>("pressure_total")),
    _Je(declareProperty<Real>("Je")),
    _s(getParam<Real>("slope")),
    _Y_final(coupledValue("Y_final")),
    _A1(getParam<Real>("A1")),
    _R1(getParam<Real>("R1")),
    _A2(getParam<Real>("A2")),
    _R2(getParam<Real>("R2")),
    _omega(getParam<Real>("omega")),
    _dMG_dT(declareProperty<Real>("dMG_dT")),
    _dJWL_dT(declareProperty<Real>("dJWL_dT")),
    _pressureHS(declareProperty<Real>("pressureHS")),
    _pressureHS2(declareProperty<Real>("pressureHS2")),
    _pressure_flag(declareProperty<Real>("pressure_flag")),
    _pressure_total_old(getMaterialPropertyOld<Real>("pressure_total")),
    _flag_threshold(getParam<Real>("flag_threshold")),
    _pressure_hist(declareProperty<Real>("pressure_hist")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _lagrangian_strain(declareProperty<RankTwoTensor>("lagrangian_strain"))
    
{
}

void
ComputeIntPVal::computeQpProperties()
{
  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour
  Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  Real ss = std::sqrt(K0 / _rho[_qp]);

  //Compute Mie Gruneisen pressure for unreacted material
  Real P_mg;
  //
  Real J = _deformation_gradient[_qp].det();
  Real J_dot = (_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det()) / _dt;
  Real eta = 1.0 - J;
  Real rho0_rho = 1.0 - eta; //express rho / rho_0
  
  RankTwoTensor epsilon_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;

  P_mg = _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J); //initial thermal expansion term
  P_mg += K0 * eta * (1.0 - (_Gamma / 2.0) * (eta)) / std::pow((1.0 - _s * eta), 2.0);
  _pressure_mg[_qp] = - P_mg; //store pressure

  //compute JWL pressure for reacted final products
  Real P_JWL = _A1 * (1.0 - _omega / (_R1 * rho0_rho)) * std::exp(- _R1 * rho0_rho); //mechanical term 1
  P_JWL += _A2 * (1.0 - _omega / (_R2 * J)) * std::exp(- _R2 * J); //mechanical term 2
  P_JWL += _omega * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) / J; //thermal expansion term
  _pressure_JWL[_qp] = - P_JWL; //store
  
  //interpolate pressure based on reaction level

  Real P_total;
  P_total = ((1.0 - _Y_final[_qp]) * (P_mg)) + (_Y_final[_qp] * P_JWL); //compute total pressure
  _pressure_total[_qp] = - P_total;

  //write total pressure into stress as a hydrostatic component
  //compute and declare the derivatives of each partial pressure wrt temperature to consume on PressureHS

  Real dMG_dT;
  Real dJWL_dT;

  dMG_dT = - _Gamma * _rho[_qp] * _Cv[_qp] * (1.0 / J);
  dJWL_dT = - _omega * _rho[_qp] * _Cv[_qp] * (1.0 / J);

  _dMG_dT[_qp] = dMG_dT;
  _dJWL_dT[_qp] = dJWL_dT;

  _pressureHS[_qp] = 0.0;
  _pressureHS2[_qp] = 0.0;

  //compute lagrangian strain and other tensors

  RankTwoTensor C = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp];
  RankTwoTensor E = 0.5 * (C - I2); //lagrangian strain
  _lagrangian_strain[_qp] = E;

  //compute pressure flag

  //compute pressure change over time

  Real dP;
  dP = (_pressure_total[_qp] - _pressure_total_old[_qp]) / _dt;
  _pressure_hist[_qp] = dP;

  if (_pressure_flag[_qp] == 0.0){ //not already calling net
    if (std::abs(_pressure_total[_qp]) >= _flag_threshold){ //pressure is increasing beyond a given threshold
      _pressure_flag[_qp] = 1.0;
    }
  }
}
