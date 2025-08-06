#include "ADComputeYdotsRDXNS.h"
#include "RankTwoTensor.h"

registerMooseObject("catApp", ADComputeYdotsRDXNS);

InputParameters
ADComputeYdotsRDXNS::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute the rate of reaction of each species for a 3 phase - 2 stage reaction model");
    params.addRequiredCoupledVar("temperature", "temperature");

    params.addRequiredParam<Real>("Z1", "pre-exponential factor Z1");
    params.addRequiredParam<Real>("Z2", "pre-exponential factor Z2");
    params.addRequiredParam<Real>("Rg", "thermodynamic constan");

    params.addRequiredParam<Real>("E1", "reaction energy 1");
    params.addRequiredParam<Real>("E2", "reaction energy 2");
    params.addRequiredParam<Real>("T_trans", "transition temperature for heat of reaction");
    params.addRequiredParam<Real>("a1", "a1 parameter for heat of reaction 1");
    params.addRequiredParam<Real>("b1", "b1 parameter for heat of reaction 1");
    params.addRequiredParam<Real>("a2", "a2 parameter for heat of reaction 2");
    params.addRequiredParam<Real>("b2", "b2 parameter for heat of reaction 2");

    params.addRequiredCoupledVar("dirac_switch_react", "dirac");
    params.addRequiredParam<Real>("switch_react", "dirac switch value to turn on reaction");
    params.addRequiredParam<Real>("rate_limit", "rate_limit");
    return params;
}

ADComputeYdotsRDXNS::ADComputeYdotsRDXNS(const InputParameters & parameters)
  : Material(parameters),
    _T(coupledValue("temperature")),

    //reaction rate parameters
    _Z1(getParam<Real>("Z1")),
    _Z2(getParam<Real>("Z2")),
    _E1(getParam<Real>("E1")),
    _E2(getParam<Real>("E2")),

    _Rg(getParam<Real>("Rg")),

    //declare rate of reaction terms
    _r1(declareADProperty<Real>("r1")),
    _r2(declareADProperty<Real>("r2")),

    _Q1(declareADProperty<Real>("Q1")),
    _Q2(declareADProperty<Real>("Q2")),
    _a1(getParam<Real>("a1")),
    _b1(getParam<Real>("b1")),
    _a2(getParam<Real>("a2")),
    _b2(getParam<Real>("b2")),
    _T_trans(getParam<Real>("T_trans")),
    _dirac_switch_react(coupledValue("dirac_switch_react")),
    _switch_react(getParam<Real>("switch_react")),
    _rate_limit(getParam<Real>("rate_limit"))
{   
}

void
ADComputeYdotsRDXNS::computeQpProperties()
{
    Real R_RDX = _Rg;
    if(_dirac_switch_react[_qp] >= _switch_react){ //turn on after reaction heat has gone by
        //compute the rates of reaction of all reactions
        _r1[_qp] = std::min(_Z1 * std::exp(- _E1 / (R_RDX * _T[_qp])), _rate_limit);
        _r2[_qp] = std::min(_Z2 * std::exp(- _E2 / (R_RDX * _T[_qp])), _rate_limit);
    }else{ //don't do anything before
        _r1[_qp] = 0.;
        _r2[_qp] = 0.;
    }

    _Q1[_qp] = (_a1) + (_b1 * std::max(0., _T[_qp] - _T_trans));
    _Q2[_qp] = (_a2) + (_b2 * std::max(0., _T[_qp] - _T_trans));  
}