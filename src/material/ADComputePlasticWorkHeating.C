#include "ADComputePlasticWorkHeating.h"
#include "RankTwoTensor.h"

registerMooseObject("catApp", ADComputePlasticWorkHeating);

InputParameters
ADComputePlasticWorkHeating::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute the plastic strain rate work term to be called by a heat source kernel");
    params.addRequiredParam<Real>("beta_p", "beta_p");
    params.addRequiredParam<bool>("use_PK2", "use_PK2");
    params.addRequiredCoupledVar("dirac_switch_react", "dirac_switch_react");
    params.addRequiredParam<Real>("thr_activation", "thr_activation");
    params.addRequiredParam<bool>("use_lump", "whether to compute lumped terms or not");
    return params;
}

ADComputePlasticWorkHeating::ADComputePlasticWorkHeating(const InputParameters & parameters)
  : Material(parameters),
    _cauchy_stress(getMaterialProperty<RankTwoTensor>("cauchy_stress")),
    _beta_p(getParam<Real>("beta_p")),
    _Ep_dot(getMaterialProperty<RankTwoTensor>("Ep_dot")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _F(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _use_PK2(getParam<bool>("use_PK2")),
    _dirac_switch_react(coupledValue("dirac_switch_react")),
    _thr_activation(getParam<Real>("thr_activation")),
    //declare heat sourves
    _q_plastic(declareADProperty<Real>("q_plastic")),
    _use_lump(getParam<bool>("use_lump"))
{   
}

void
ADComputePlasticWorkHeating::computeQpProperties()
{
    RankTwoTensor I2(RankTwoTensor::initIdentity);
    ADReal q_plastic;

    //whether to use PK2 or cauchy stress
    if (_use_PK2){
        Real J = _F[_qp].det();
        RankTwoTensor PK2 = J * _F[_qp].inverse() * _cauchy_stress[_qp] * _F[_qp].inverse().transpose();
        q_plastic = std::max(PK2.doubleContraction(_Ep_dot[_qp]), 0.);
    }else{
        q_plastic = _cauchy_stress[_qp].doubleContraction(_Ep_dot[_qp]);
    }

    //always multiply by beta_p
    q_plastic *= _beta_p;

    //whether to use lump constant or not
    if (_use_lump){
        q_plastic *= (1. / (_rho[_qp] * _cv[_qp]));
    }

    //activation for MISTERnet simulation
    if(_dirac_switch_react[_qp] > _thr_activation){
        q_plastic *= 1.; //unchanged, just forwards and outputs
    }else{
        q_plastic *= 0.; //if the prediction isn't over, the heat from plasticity is zero
    }

    //append into actual material property
    _q_plastic[_qp] = q_plastic;
}