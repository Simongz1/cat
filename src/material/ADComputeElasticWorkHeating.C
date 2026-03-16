#include "ADComputeElasticWorkHeating.h"

registerMooseObject("mlApp", ADComputeElasticWorkHeating);

//test: currently changing T to retrieve ADVariableValue

InputParameters
ADComputeElasticWorkHeating::validParams()
{
    InputParameters params = DerivativeMaterialInterface<ADMaterial>::validParams();
    params.addClassDescription("compute the elastic compression work term to be called by a heat source kernel");
    params.addCoupledVar("temperature", "temperature");
    params.addRequiredParam<Real>("beta_av", "beta_av");
    params.addRequiredCoupledVar("dirac_switch_react", "dirac_switch_react");
    params.addRequiredParam<Real>("thr_activation", "thr_activation");
    params.addParam<bool>("consistent_tau", true, "consistent_tau");
    return params;
}

ADComputeElasticWorkHeating::ADComputeElasticWorkHeating(const InputParameters & parameters)
  : DerivativeMaterialInterface<ADMaterial>(parameters),
    _T(adCoupledValue("temperature")),
    _Tgrad(adCoupledGradient("temperature")),
    _beta_av(getParam<Real>("beta_av")),
    _S(getADMaterialProperty<RankTwoTensor>("S")), //this is PK2 Stress
    _dPdT(getADMaterialProperty<Real>("dPdT")),
    _Ee_dot(getADMaterialProperty<RankTwoTensor>("Ee_dot")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _dirac_switch_react(adCoupledValue("dirac_switch_react")),
    _Fe(getADMaterialProperty<RankTwoTensor>("Fe")),
    _thr_activation(getParam<Real>("thr_activation")),
    _q_elastic(declareADProperty<Real>("q_elastic")),
    _norm_gradT(declareADProperty<Real>("norm_gradT")),
    _consistent_tau(getParam<bool>("consistent_tau")),
    _time_react(getADMaterialProperty<Real>("time_react"))
{   
}

void
ADComputeElasticWorkHeating::computeQpProperties()
{
    ADRankTwoTensor I2;
    I2.setToIdentity();

    //component contribution from volumetric compression
    ADReal q_pressure;
    ADRankTwoTensor Ce = _Fe[_qp].transpose() * _Fe[_qp];

    //if use PK2, use work conjugate C.inverse()
    q_pressure = - std::max(_T[_qp] * _dPdT[_qp] * (Ce.inverse().doubleContraction(_Ee_dot[_qp])), ADReal(0.0));

    ADReal reaction_thr;
    reaction_thr = _consistent_tau ? _time_react[_qp] : ADReal(_thr_activation);

    //activation for MISTERnet simulations
    if(_dirac_switch_react[_qp] > reaction_thr){
        q_pressure *= 1.; //keep while activated
    }else{
        q_pressure *= 0.; //set to zero before activation
    }

    _q_elastic[_qp] = q_pressure;
    _norm_gradT[_qp] = _Tgrad[_qp].norm();
}