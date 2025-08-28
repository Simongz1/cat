#include "ADComputeElasticWorkHeating.h"

registerMooseObject("catApp", ADComputeElasticWorkHeating);

//test: currently changing T to retrieve ADVariableValue

InputParameters
ADComputeElasticWorkHeating::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute the elastic compression work term to be called by a heat source kernel");
    params.addCoupledVar("temperature", "temperature");
    params.addRequiredParam<Real>("beta_av", "beta_av");
    params.addRequiredCoupledVar("dirac_switch_react", "dirac_switch_react");
    params.addRequiredParam<Real>("thr_activation", "thr_activation");
    params.addRequiredParam<bool>("use_PK2", "whether to use PK2 work conjugate configuration");
    params.addRequiredParam<bool>("use_lump", "whether to compute lumped terms or not");
    return params;
}

ADComputeElasticWorkHeating::ADComputeElasticWorkHeating(const InputParameters & parameters)
  : Material(parameters),
    _T(adCoupledValue("temperature")),
    _beta_av(getParam<Real>("beta_av")),
    _pressure_av(getADMaterialProperty<Real>("pressure_av")),
    _dP_dT(getADMaterialProperty<Real>("dP_dT")),
    _Ee_dot(getMaterialProperty<RankTwoTensor>("Ee_dot")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _dirac_switch_react(coupledValue("dirac_switch_react")),
    _Fe(getMaterialProperty<RankTwoTensor>("Fe")),
    _thr_activation(getParam<Real>("thr_activation")),
    _use_PK2(getParam<bool>("use_PK2")),
    //declare properties
    _q_elastic(declareADProperty<Real>("q_elastic")),
    _use_lump(getParam<bool>("use_lump"))
{   
}

void
ADComputeElasticWorkHeating::computeQpProperties()
{
    RankTwoTensor I2(RankTwoTensor::initIdentity);

    //component contribution from volumetric compression
    ADReal q_pressure;
    RankTwoTensor Ce = _Fe[_qp].transpose() * _Fe[_qp];

    //if use PK2, use work conjugate C.inverse()
    if (_use_PK2){
        q_pressure = - std::max(_T[_qp] * _dP_dT[_qp] * (Ce.inverse().doubleContraction(_Ee_dot[_qp])), 0.);
    }else{
        q_pressure = - std::max(_T[_qp] * _dP_dT[_qp] * _Ee_dot[_qp].trace(), 0.0);
    }

    //add artificial viscosity factor

    ADReal q_av;
    q_av = _beta_av * _pressure_av[_qp] * _Ee_dot[_qp].trace();

    ADReal q_tot = q_pressure + q_av;

    //activation for MISTERnet simulations
    if(_dirac_switch_react[_qp] > _thr_activation){
        q_tot *= 1.; //keep while activated
    }else{
        q_tot *= 0.; //set to zero before activation
    }

    _q_elastic[_qp] = q_tot;
}