#include "ComputeIntPValRDX.h"

registerMooseObject("TensorMechanicsApp", ComputeIntPValRDX);

InputParameters
ComputeIntPValRDX::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Standard compute Mie Gruneisen Pressure with JWL pressure for reacted material. Also computes artificial viscosity contribution");
  params.addRequiredParam<Real>("T_ref", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("element_size", "element_size");
  params.addRequiredCoupledVar("Y_final", "final products mass fraction");
  params.addRequiredParam<Real>("A_u", "JWL reacted EOS parameter A1");
  params.addRequiredParam<Real>("R1_u", "JWL reacted EOS parameter B1");
  params.addRequiredParam<Real>("B_u", "JWL reacted EOS parameter A2");
  params.addRequiredParam<Real>("R2_u", "JWL reacted EOS parameter B2");
  params.addRequiredParam<Real>("A_r", "JWL reacted EOS parameter A1");
  params.addRequiredParam<Real>("R1_r", "JWL reacted EOS parameter B1");
  params.addRequiredParam<Real>("B_r", "JWL reacted EOS parameter A2");
  params.addRequiredParam<Real>("R2_r", "JWL reacted EOS parameter B2");
  params.addRequiredParam<Real>("omega_u", "JWL reacted EOS parameter omega");
  params.addRequiredParam<Real>("omega_r", "JWL reacted EOS parameter omega");
  params.addRequiredParam<Real>("flag_threshold", "pressure flag threshold value");
  params.addRequiredParam<Real>("P0", "pressure at ambient temperature");
  params.addRequiredParam<Real>("use_RDX", "use RDX model or HMX");
  //request HMX parameters
  params.addRequiredParam<Real>("Gamma", "Grunseisen Parameter");
  params.addRequiredParam<Real>("slope", "slope");
  params.addRequiredParam<Real>("A1", "JWL reacted EOS parameter A1");
  params.addRequiredParam<Real>("R1", "JWL reacted EOS parameter B1");
  params.addRequiredParam<Real>("A2", "JWL reacted EOS parameter A2");
  params.addRequiredParam<Real>("R2", "JWL reacted EOS parameter B2");
  params.addRequiredParam<Real>("omega", "JWL reacted EOS parameter omega");
  //velocity for calling mistnet
  params.addRequiredCoupledVar("vx", "x component of velocity");
  params.addRequiredCoupledVar("ax", "x component of acceleration");
  params.addRequiredParam<Real>("thr_a", "acceleration threshold");
  params.addRequiredParam<Real>("thr_v", "velocity threshold");
  return params;
}

ComputeIntPValRDX::ComputeIntPValRDX(
    const InputParameters & parameters)
  : Material(parameters),
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
    _Y_final(coupledValue("Y_final")),
    _A_u(getParam<Real>("A_u")),
    _R1_u(getParam<Real>("R1_u")),
    _B_u(getParam<Real>("B_u")),
    _R2_u(getParam<Real>("R2_u")),
    _omega_u(getParam<Real>("omega_u")),
    _A_r(getParam<Real>("A_r")),
    _R1_r(getParam<Real>("R1_r")),
    _B_r(getParam<Real>("B_r")),
    _R2_r(getParam<Real>("R2_r")),
    _omega_r(getParam<Real>("omega_r")),
    _dMG_dT(declareProperty<RankTwoTensor>("dMG_dT")),
    _dJWL_dT(declareProperty<RankTwoTensor>("dJWL_dT")),
    _pressure_total_old(getMaterialPropertyOld<Real>("pressure_total")),
    _pressure_hist(declareProperty<Real>("pressure_hist")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _f_inv(getMaterialProperty<RankTwoTensor>("inverse_incremental_deformation_gradient")),
    _P0(getParam<Real>("P0")),
    _use_RDX(getParam<Real>("use_RDX")),
    //get HMX parameters
    _Gamma(getParam<Real>("Gamma")),
    _s(getParam<Real>("slope")),
    _A1(getParam<Real>("A1")),
    _R1(getParam<Real>("R1")),
    _A2(getParam<Real>("A2")),
    _R2(getParam<Real>("R2")),
    _omega(getParam<Real>("omega")),
    _vx(coupledValue("vx")),
    _ax(coupledValue("ax")),
    _thr_a(getParam<Real>("thr_a")),
    _thr_v(getParam<Real>("thr_v")),
    _v_flag(declareProperty<Real>("v_flag")),
    _J(declareProperty<Real>("J"))
{
}

void 
ComputeIntPValRDX::initQpStatefulProperties()
{
  _v_flag[_qp] = 0.0; //initialize flag
}
void
ComputeIntPValRDX::computeQpProperties()
{
  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour
  Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  Real ss = std::sqrt(K0 / _rho[_qp]);

  //Compute Mie Gruneisen pressure for unreacted material
  Real P_mg;
  Real P_JWL;
  //
  Real J = _deformation_gradient[_qp].det();
  _J[_qp] = J;
  Real J_dot = (_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det()) / _dt;
  Real eta = 1.0 - J;
  Real rho0_rho = 1.0 - eta; //express rho / rho_0
  
  RankTwoTensor epsilon_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;

  P_mg = _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J); //initial thermal expansion term
  P_mg += K0 * eta * (1.0 - (_Gamma / 2.0) * (eta)) / std::pow((1.0 - _s * eta), 2.0);
  _pressure_mg[_qp] = - P_mg; //store pressure 

  P_JWL = _A1 * std::exp(- _R1 * J) + _A2 * std::exp(- _R2 * J) + _omega * _rho[_qp] * _Cv[_qp] * _T[_qp] / J;
  _pressure_JWL[_qp] = - P_JWL;

  //P_JWL = _A1 * (1.0 - _omega / (_R1 * J)) * std::exp(- _R1 * J); //mechanical term 1
  //P_JWL += _A2 * (1.0 - _omega / (_R2 * J)) * std::exp(- _R2 * J); //mechanical term 2
  //P_JWL += _omega * _rho[_qp] * _Cv[_qp] * (_T[_qp]) / J; //thermal expansion term
  //_pressure_JWL[_qp] = - P_JWL; //store

  //P_mg = _A_u * (1.0 - (_omega_u / (_R1_u * J))) * std::exp(- _R1_u * (J));
  //P_mg += _B_u * (1.0 - (_omega_u / (_R2_u * J))) * std::exp(- _R2_u * (J));
  //P_mg += _omega_u * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J);
  //P_mg += _P0;
  //_pressure_mg[_qp] = - P_mg; //store pressure

  //compute JWL pressure for reacted final products
  //P_JWL = _A_r * (1.0 - (_omega_r / (_R1_r * J))) * std::exp(- _R1_r * (J));
  //P_JWL += _B_r * (1.0 - (_omega_r / (_R2_r * J))) * std::exp(- _R2_r * (J));
  //P_JWL += _omega_r * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J);
  //_pressure_JWL[_qp] = - P_JWL; //store

  //if (_use_RDX == 1.0){
  //  Real P_mg = _A_u * (1.0 - (_omega_u / (_R1_u * J))) * std::exp(- _R1_u * (J));
  //  P_mg += _B_u * (1.0 - (_omega_u / (_R2_u * J))) * std::exp(- _R2_u * (J));
  //  P_mg += _omega_u * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J);
  //  P_mg += _P0;
  //  _pressure_mg[_qp] = - P_mg; //store pressure
//
  //  //compute JWL pressure for reacted final products
  //  Real P_JWL = _A_r * (1.0 - (_omega_r / (_R1_r * J))) * std::exp(- _R1_r * (J));
  //  P_JWL += _B_r * (1.0 - (_omega_r / (_R2_r * J))) * std::exp(- _R2_r * (J));
  //  P_JWL += _omega_r * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J);
  //  _pressure_JWL[_qp] = - P_JWL; //store
  //}

  //interpolate pressure based on reaction level
  //if (_use_RDX == 0.0){
  //  //UNREACTED PRESSURE FOR HMX
  //  Real P_mg = _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J); //initial thermal expansion term
  //  //P_mg += K0 * eta * (1.0 - (_Gamma / 2.0) * (eta)) / std::pow((1.0 - _s * eta), 2.0);
  //  P_mg += ((K0 * eta) / std::pow((1.0 - _s * eta), 2.0)) * ((_Gamma / 2.0) * (J - 1) - 1);
  //  _pressure_mg[_qp] = - P_mg; //store pressure
//
  //  //REACTED PRESSURE FOR HMX
  //  Real P_JWL = _A1 * std::exp(- _R1 * J); //mechanical term 1
  //  P_JWL += _A2 * std::exp(- _R2 * J); //mechanical term 2
  //  P_JWL += _omega * _rho[_qp] * _Cv[_qp] * (_T[_qp]); //thermal expansion term
  //  _pressure_JWL[_qp] = - P_JWL; //store
  //}
  
  //pressure interpolation
  Real P_total;
  P_total = ((1.0 - _Y_final[_qp]) * (P_mg)) + (_Y_final[_qp] * P_JWL); //compute total pressure
  _pressure_total[_qp] = - P_total;

  //write total pressure into stress as a hydrostatic component
  //compute and declare the derivatives of each partial pressure wrt temperature to consume on PressureHS

  //compression heat source declaration
  RankTwoTensor dMG_dT;
  RankTwoTensor dJWL_dT;

  dMG_dT = - _omega_u * _rho[_qp] * _Cv[_qp] * (1.0 / J);
  dJWL_dT = - _omega_r * _rho[_qp] * _Cv[_qp] * (1.0 / J);

  _dMG_dT[_qp] = dMG_dT;
  _dJWL_dT[_qp] = dJWL_dT;

  //compute lagrangian strain and other tensors
  RankTwoTensor f = _f_inv[_qp].inverse();
  RankTwoTensor f_bar = f / std::cbrt(f.det());
  
  //conditions: (a < value) and (v > value)
  //v_flag == 0 means that the flag has not been set to 1 here
  if (_v_flag[_qp] == 0.0 && _ax[_qp] <= _thr_a && _vx[_qp] >= _thr_v){
    _v_flag[_qp] = 1.0; //set flag to one, call in misternet material as condition for ComputeQpProperties
  }
}
