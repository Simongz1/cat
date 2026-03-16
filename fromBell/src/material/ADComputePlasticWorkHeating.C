#include "ADComputePlasticWorkHeating.h"
#include "RankTwoTensor.h"

registerMooseObject("mistApp", ADComputePlasticWorkHeating);

InputParameters
ADComputePlasticWorkHeating::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute the plastic strain rate work term to be called by a heat source kernel");
    params.addRequiredParam<Real>("beta_p", "beta_p");
    params.addRequiredCoupledVar("dirac_switch_react", "dirac_switch_react");
    params.addRequiredParam<Real>("thr_activation", "thr_activation");
    return params;
}

ADComputePlasticWorkHeating::ADComputePlasticWorkHeating(const InputParameters & parameters)
  : Material(parameters),
    _S(getADMaterialProperty<RankTwoTensor>("S")),
    _beta_p(getParam<Real>("beta_p")),
    _Ep_dot(getADMaterialProperty<RankTwoTensor>("Ep_dot")),
    _Je(getADMaterialProperty<Real>("Je")),
    _dirac_switch_react(adCoupledValue("dirac_switch_react")),
    _thr_activation(getParam<Real>("thr_activation")),
    //declare heat sourves
    _q_plastic(declareADProperty<Real>("q_plastic"))
{   
}

void
ADComputePlasticWorkHeating::computeQpProperties()
{
    RankTwoTensor I2(RankTwoTensor::initIdentity);
    ADReal q_plastic;

    q_plastic = _beta_p * std::max(_S[_qp].doubleContraction(_Ep_dot[_qp]), ADReal(0.0));

    //activation for MISTERnet simulation
    if(_dirac_switch_react[_qp] > _thr_activation){
        q_plastic *= 1.; //unchanged, just forwards and outputs
    }else{
        q_plastic *= 0.; //if the prediction isn't over, the heat from plasticity is zero
    }

    //append into actual material property
    _q_plastic[_qp] = q_plastic;
}