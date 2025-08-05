#include "ComputeInterpolatedPressure.h"

registerMooseObject("TensorMechanicsApp", ComputeInterpolatedPressure);

InputParameters
ComputeInterpolatedPressure::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
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

ComputeInterpolatedPressure::ComputeInterpolatedPressure(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
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
    _pressure_av(declareProperty<Real>("pressure_av")),
    _pressureHS(declareProperty<Real>("pressureHS")),
    _pressureHS2(declareProperty<Real>("pressureHS2")),
    _pressure_flag(declareProperty<Real>("pressure_flag")),
    _pressure_total_old(getMaterialPropertyOld<Real>("pressure_total")),
    _flag_threshold(getParam<Real>("flag_threshold")),
    _pressure_hist(declareProperty<Real>("pressure_hist"))
    
{
}

void
ComputeInterpolatedPressure::initQpStatefulProperties()
{
	_stress[_qp].zero();
  _pressure_flag[_qp] = 0.0; //initialize pressure flag
  _pressure_hist[_qp] = 0.0; //initialize pressure history variable
}

void
ComputeInterpolatedPressure::computeQpStress()
{
  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour
  Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  Real ss = std::sqrt(K0 / _rho[_qp]);

  //Compute Mie Gruneisen pressure for unreacted material
  Real P_mg;
  Real eta = - _mechanical_strain[_qp].trace();
  Real rho0_rho = 1.0 - eta; //express rho / rho_0
  
  RankTwoTensor epsilon_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  Real J_dot = epsilon_dot.trace();

  P_mg = _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / (1.0 - eta)); //initial thermal expansion term
  P_mg += K0 * eta * (1.0 - (_Gamma / 2.0) * (1.0 / rho0_rho)) / std::pow((1.0 - _s * eta), 2.0);
  _pressure_mg[_qp] = - P_mg; //store pressure

  //compute JWL pressure for reacted final products
  Real P_JWL = _A1 * (1.0 - _omega / (_R1 * rho0_rho)) * std::exp(- _R1 * rho0_rho); //mechanical term 1
  P_JWL += _A2 * (1.0 - _omega / (_R2 * rho0_rho)) * std::exp(- _R2 * rho0_rho); //mechanical term 2
  P_JWL += _omega * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) / rho0_rho; //thermal expansion term
  _pressure_JWL[_qp] = - P_JWL; //store
  
  //interpolate pressure based on reaction level

  Real P_total;
  P_total = ((1.0 - _Y_final[_qp]) * (P_mg)) + (_Y_final[_qp] * P_JWL); //compute total pressure
  _pressure_total[_qp] = - P_total;

  //write total pressure into stress as a hydrostatic component
  _stress[_qp] = - P_total * I2;
	
  //Compute artificial viscosity term
  Real P_av;
  Real Je_dot = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;

  P_av = _C0 * _rho[_qp] * (rho0_rho * std::abs(Je_dot) / std::pow(rho0_rho, 2.0)) * std::pow(_Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (Je_dot / rho0_rho) * _Le;
  _pressure_av[_qp] = P_av;
  _stress[_qp] += P_av * I2;

  //compute and declare the derivatives of each partial pressure wrt temperature to consume on PressureHS

  Real dMG_dT;
  Real dJWL_dT;

  dMG_dT = - _Gamma * _rho[_qp] * _Cv[_qp] * (1.0 / (1.0 - eta));
  dJWL_dT = - _omega * _rho[_qp] * _Cv[_qp] * 1.0 / rho0_rho;

  _dMG_dT[_qp] = dMG_dT;
  _dJWL_dT[_qp] = dJWL_dT;
  RankTwoTensor epsilon_e_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  Real epsilon_e_trace_dot = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;

  _pressureHS[_qp] = _dMG_dT[_qp] * _T[_qp] * epsilon_e_dot.trace();
  _pressureHS2[_qp] = _dMG_dT[_qp] * _T[_qp] * epsilon_e_trace_dot;

  //compute pressure flag

  //compute pressure change over time

  Real dP;
  dP = (_pressure_total[_qp] - _pressure_total_old[_qp]) / _dt;
  _pressure_hist[_qp] = dP;

  //if (_pressure_flag[_qp] == 0.0){
  //  if (dP >= _flag_threshold && std::abs(_pressure_total[_qp]) <= std::abs(_pressure_total_old[_qp])){
  //    _pressure_flag[_qp] = 0.0; //case 1. pressure still increasing
  //  }

  //  if (dP >= _flag_threshold && std::abs(_pressure_total[_qp]) >= std::abs(_pressure_total_old[_qp])){
  //    _pressure_flag[_qp] = 1.0; //case 2. maximum pressure achieved, call misternet
  //  }
  //}

  if (_pressure_flag[_qp] == 0.0){ //not already calling net
    if (std::abs(_pressure_total[_qp]) >= _flag_threshold){ //pressure is increasing beyond a given threshold
      _pressure_flag[_qp] = 1.0;
    }
  }
}
