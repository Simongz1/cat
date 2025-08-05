#include "VisHS.h"

registerMooseObject("catApp", VisHS);

InputParameters
VisHS::validParams()
{
  InputParameters params = HeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute heating from viscous dissipation");
  params.addRequiredCoupledVar("temperature", "temperature variable");
  params.addRequiredCoupledVar("c", "c");
  params.addRequiredParam<Real>("beta_av", "beta_av");
  params.addRequiredParam<Real>("beta_p", "beta_p");
  params.addRequiredParam<Real>("shock_heat", "flag to use shock heat for lipit");
  return params;
}

VisHS::VisHS(const InputParameters & parameters)
  : HeatSource(parameters),
    //retrieve Ydots
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _cauchy_stress(getMaterialProperty<RankTwoTensor>("cauchy_stress")),
    _density(getMaterialProperty<Real>("density")),
    _pressure_av(getMaterialProperty<RankTwoTensor>("pressure_av")),
    _be(getMaterialProperty<RankTwoTensor>("volume_preserving_elastic_left_cauchy_green_strain")), //plastic left stretch tensor
    _be_old(getMaterialPropertyOld<RankTwoTensor>("volume_preserving_elastic_left_cauchy_green_strain")), //old plastic left stretch tensor
    _beta_av(getParam<Real>("beta_av")),
    _beta_p(getParam<Real>("beta_p")),
    _lagrangian_strain_rate(getMaterialProperty<RankTwoTensor>("lagrangian_strain_rate")),
    _ep_rate(getMaterialProperty<Real>("ep_rate")),
    _Np(getMaterialProperty<RankTwoTensor>("flow_direction")),
    _plastic_strain(getMaterialProperty<RankTwoTensor>("plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("plastic_strain")),
    _Cijkl(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _alpha(getMaterialProperty<Real>("alpha")),
    _temperature(coupledValue("temperature")),
    _shock_heat(getParam<Real>("shock_heat")),
    _c(coupledValue("c"))
{}

Real
VisHS::computeQpResidual()
{
  RankTwoTensor F, C, E, L, d, F_dot, E_dot;
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  F = _deformation_gradient[_qp];
  F_dot = (_deformation_gradient[_qp] - _deformation_gradient_old[_qp]) / _dt;
  C = F.transpose() * F;
  E = (1.0 / 2.0) * (C - I2);
  E_dot = (1.0 / 2.0) * ((F_dot.transpose() * F_dot) - I2);
  L = F_dot * F.inverse();
  d = (1.0 / 2.0) * (L + L.transpose());
  Real J = F.det();

  //compute viscous part of stress
  RankTwoTensor sigma_vis;
  sigma_vis = _cauchy_stress[_qp] - (1.0 / 3.0) * _cauchy_stress[_qp].trace() * I2;
  RankTwoTensor be_rate = (_be[_qp] - _be_old[_qp]) / _dt;

  
  //only artificial viscosity contribution is considered
  

  //get elastic-plastic decomposition from converged Be'
  RankTwoTensor Be = std::pow(J, 4. / 3.) * _be[_qp];
  RankTwoTensor Ce = Be.transpose();
  RankTwoTensor Ee = 0.5 * (Ce - I2);

  //old
  Real J_old = _deformation_gradient_old[_qp].det();
  RankTwoTensor Be_old = std::pow(J_old, 4. / 3.) * _be_old[_qp];
  RankTwoTensor Ce_old = Be_old.transpose();
  RankTwoTensor Ee_old = 0.5 * (Ce_old - I2);

  //elastic strain
  RankTwoTensor Ee_rate = (Ee - Ee_old) / _dt;
  RankTwoTensor Ep = E - Ee;
  RankTwoTensor Ep_dot = E_dot - Ee_rate;

  //write expression
  Real q_vis_av = 0.0;
  RankTwoTensor Ep_rate;
  Ep_rate = _ep_rate[_qp] * I2; //testing Ep_dot = ep * flowdirection
  q_vis_av = _beta_av * _pressure_av[_qp].doubleContraction(_lagrangian_strain_rate[_qp]);

  RankTwoTensor cauchy_pos;
  RankTwoTensor cauchy_copy = _cauchy_stress[_qp];
  
  unsigned int i, j;

  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      if (cauchy_copy(i, j) < 0.){
        cauchy_pos(i, j) = 0.;
      }
      else{
        cauchy_pos(i, j) = cauchy_copy(i, j);
      }
    }
  }

  Real q_vis = _beta_p * cauchy_pos.doubleContraction(Ep_rate); //NOTE: originally plastic work was 
  RankTwoTensor prate = (_plastic_strain[_qp] - _plastic_strain_old[_qp]) / _dt;
  Real q_vis_p = std::abs(_beta_p * _cauchy_stress[_qp].doubleContraction(prate));
  Real q_vis_fromep = _beta_p * _cauchy_stress[_qp].doubleContraction(Ep_rate);

  //include compression heat source, specifically for LIPIT

  const Real K = ElasticityTensorTools::getIsotropicBulkModulus(_Cijkl[_qp]);

  RankTwoTensor Cinv = Ce.inverse();
  Real q_comp = - _alpha[_qp] * K * _temperature[_qp] * (Cinv.doubleContraction(Ee_rate));
  //return - q_vis_p * _test[_i][_qp];

  return - q_vis_p * _test[_i][_qp];
}

Real
VisHS::computeQpJacobian()
{
  Real Jac = 0.0;
  return Jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
VisHS::computeQpOffDiagonalJacobian(unsigned int jvar)
{
  Real OffJac;
  OffJac = 0.0;
  return OffJac * _test[_i][_qp] * _phi[_j][_qp];
}
