#include "PressureHS.h"

registerMooseObject("catApp", PressureHS);

InputParameters
PressureHS::validParams()
{
  InputParameters params = HeatSource::validParams();
  params.addClassDescription("compute elastic thermomechanical coupling heat source term");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredCoupledVar("Yfinal", "reacted mass fraction");
  params.addRequiredParam<Real>("beta_av", "beta_av");
  return params;
}

PressureHS::PressureHS(const InputParameters & parameters)
  : HeatSource(parameters),
    _T(coupledValue("temperature")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>("mechanical_strain")),
    _Yfinal(coupledValue("Yfinal")),
    _dMG_dT(getMaterialProperty<RankTwoTensor>("dMG_dT")),
    _dJWL_dT(getMaterialProperty<RankTwoTensor>("dJWL_dT")),
    _beta_av(getParam<Real>("beta_av")),
    _pressure_av(getMaterialProperty<RankTwoTensor>("pressure_av")),
    _lagrangian_strain_rate(getMaterialProperty<RankTwoTensor>("lagrangian_strain"))
{}

Real
PressureHS::computeQpResidual()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  
  //compute MG contribution
  RankTwoTensor epsilon_e_dot;
  Real epsilon_e_trace_dot;
  epsilon_e_trace_dot = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  epsilon_e_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  //
  Real res_MG;
  res_MG = _u[_qp] * _dMG_dT[_qp].doubleContraction(_lagrangian_strain_rate[_qp]); //contribution of MG to heat term

  //contrubtion of JWL
  Real res_JWL;
  res_JWL = _u[_qp] * _dJWL_dT[_qp].doubleContraction(_lagrangian_strain_rate[_qp]);

  //add as an interpolation for both reacted and unreacted EOS

  Real res_pressure;
  res_pressure = ((1.0 - _Yfinal[_qp]) * (res_MG)) + (_Yfinal[_qp] * res_JWL);

  //Compute artificial viscosity term
  
  Real res_av;
  RankTwoTensor P_av;
  P_av = _pressure_av[_qp];
  res_av = _beta_av * P_av.doubleContraction(epsilon_e_dot);
  Real res_tot = res_pressure + res_av;
  return - res_tot * _test[_i][_qp];

}

Real
PressureHS::computeQpJacobian()
{
  Real Jac_MG;
  Real Jac_JWL;
  Real Jac_av;
  Real Jac_tot;

  Jac_MG = _dMG_dT[_qp].trace();
  Jac_av = 0.0;
  Jac_JWL = _dJWL_dT[_qp].trace();
  Jac_tot = ((1. - _Yfinal[_qp]) * Jac_MG) + (_Yfinal[_qp] * Jac_JWL) + Jac_av;
  return - Jac_tot * _test[_i][_qp] * _phi[_j][_qp];
}

Real
PressureHS::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  Real OffJac;
  OffJac = 0.0;
  return OffJac * _test[_i][_qp] * _phi[_j][_qp];
}