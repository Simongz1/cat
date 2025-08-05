#include "MieGruneisenStress.h"

registerMooseObject("TensorMechanicsApp", MieGruneisenStress);

InputParameters
MieGruneisenStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Compute Elastic MieGruneisen EOS Pressure with Thermomechanical Coupling");
  params.addRequiredParam<Real>("Gamma", "Grunseisen Parameter");
  params.addRequiredParam<Real>("T_ref", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("viscosity_type", "viscosity_type");
  params.addRequiredParam<Real>("element_size", "element_size");
  params.addRequiredParam<Real>("slope", "slope");
  params.addRequiredParam<Real>("A", "A");
  params.addRequiredParam<Real>("B", "B");
  params.addRequiredParam<Real>("R1", "R1");
  params.addRequiredParam<Real>("R2", "R2");
  params.addRequiredParam<Real>("omega", "omega");
  params.addRequiredCoupledVar("Y4", "Y4");
  params.addRequiredParam<Real>("use_JWL", "use_JWL");
  return params;
}

MieGruneisenStress::MieGruneisenStress(
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
    _current_elem(_assembly.elem()),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _viscosity_type(getParam<Real>("viscosity_type")),
    _Le(getParam<Real>("element_size")),
    _pressure_mg(declareProperty<Real>("pressure_mg")),
    _pressure_JWL(declareProperty<Real>("pressure_JWL")),
    _pressure_total(declareProperty<Real>("pressure_total")),
    _pressure_av(declareProperty<Real>("pressure_av")),
    _s(getParam<Real>("slope")),
    _A(getParam<Real>("A")),
    _B(getParam<Real>("B")),
    _R1(getParam<Real>("R1")),
    _R2(getParam<Real>("R2")),
    _omega(getParam<Real>("omega")),
    _Y4(coupledValue("Y4")),
    _use_JWL(getParam<Real>("use_JWL")),
    _lambda(declareProperty<Real>("lambda")),
    _rhorho(declareProperty<Real>("rhorho")),
    _dMG_dT(declareProperty<Real>("dMG_dT")),
    _dJWL_dT(declareProperty<Real>("dJWL_dT"))
{
}

void
MieGruneisenStress::initQpStatefulProperties()
{
	_stress[_qp].zero();
}

void
MieGruneisenStress::computeQpStress()
{
  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  Real ss = std::sqrt(K0 / _rho[_qp]);

  //Compute Mie Gruneisen pressure for unreacted material
  Real P_mg;
  Real eta = - _mechanical_strain[_qp].trace();
  Real rho0_rho = 1.0 - eta; //express rho / rho_0
  _rhorho[_qp] = rho0_rho;
  
  RankTwoTensor epsilon_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  Real J_dot = epsilon_dot.trace();

  P_mg = _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / (1.0 - eta)); //initial thermal expansion term
  P_mg += K0 * eta * (1.0 - (_Gamma / 2.0) * (1.0 / rho0_rho)) / std::pow((1.0 - _s * eta), 2.0);
  _pressure_mg[_qp] = - P_mg; //store pressure
  
  ///compute JWL pressure for reacted final products
  Real P_JWL = _A * (1.0 - _omega / (_R1 * rho0_rho)) * std::exp(- _R1 * rho0_rho); //mechanical term 1
  P_JWL += _B * (1.0 - _omega / (_R2 * rho0_rho)) * std::exp(- _R2 * rho0_rho); //mechanical term 2
  P_JWL += _omega * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) / rho0_rho; //thermal expansion term
  _pressure_JWL[_qp] = - P_JWL; //store in material property
  
  //declare lambda for use in multiple materials / kernels
  
  Real P_total;
  
  //write partial pressures into stress tensor
  if (_use_JWL == 1)
  {
  	P_total = ((1.0 - _Y4[_qp]) * P_mg) + (_Y4[_qp] * P_JWL);

  }
  else
  {
  	P_total = P_mg;
  }
  
  _pressure_total[_qp] = - P_total;
  
  //write into stress tensor as volumetric contribution
  _stress[_qp] = - P_total * I2;
	
  //Compute artificial viscosity term
  Real P_av;
  Real Le;
  Real Je_dot = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;

  //assign element size for viscosity
  if(_viscosity_type == 1){ //representing constant
    Le = _Le;
  }
  else {
    Le = _current_elem->hmax();
  }
  
  P_av = _C0 * _rho[_qp] * (rho0_rho * std::abs(Je_dot) / std::pow(rho0_rho, 2.0)) * std::pow(_Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (Je_dot / rho0_rho) * _Le;
  _pressure_av[_qp] = P_av;
  
  _stress[_qp] += P_av * I2; //edited sound speed CHECK
  
  //compute and declare the derivatives of each partial pressure wrt temperature to consume on PressureHS

  Real dMG_dT;
  Real dJWL_dT;

  dMG_dT = _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp]) * (1.0 / (1.0 - eta));
  dJWL_dT = _omega * _rho[_qp] * _Cv[_qp] * (_T[_qp]) / rho0_rho;

  _dMG_dT[_qp] = dMG_dT;
  _dJWL_dT[_qp] = dJWL_dT;
}
