#include "ADComputeYdotsRDXNS.h"
#include "RankTwoTensor.h"

registerMooseObject("catApp", ADComputeYdotsRDXNS);

//test: currently changing T to retrieve ADVariableValue

InputParameters
ADComputeYdotsRDXNS::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute the rate of reaction of each species for a 3 phase - 2 stage reaction model");
    params.addCoupledVar("temperature", "temperature");

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

    //parameters to compute decomposition heat rate

    params.addCoupledVar("Y1", "Y1");
    params.addCoupledVar("Y2", "Y2");
    params.addRequiredParam<bool>("use_lump", "whether to compute lumped terms or not");
    params.addRequiredParam<bool>("dynamic_tau", "use dynamic tau");
    params.addRequiredParam<Real>("thr_activation_rates", "thr_activation_rates");
    return params;
}

ADComputeYdotsRDXNS::ADComputeYdotsRDXNS(const InputParameters & parameters)
  : Material(parameters),
    _T(adCoupledValue("temperature")),

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
    _Y1(adCoupledValue("Y1")),
    _Y2(adCoupledValue("Y2")),

    _use_lump(getParam<bool>("use_lump")),
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    //declarations for kernels
    _q_decomposition(declareADProperty<Real>("q_decomposition")),
    _Y1_dot(declareADProperty<Real>("Y1_dot")),
    _Y2_dot(declareADProperty<Real>("Y2_dot")),
    _Y3_dot(declareADProperty<Real>("Y3_dot")),
    _dynamic_tau(getParam<bool>("dynamic_tau")),
    _time_react(getMaterialProperty<Real>("time_react")),
    _thr(getParam<Real>("thr_activation_rates"))
{   
}

void
ADComputeYdotsRDXNS::computeQpProperties()
{
    Real R_RDX = _Rg;

    bool condition_chemistry;
    Real cutoff;
    if (_dynamic_tau){
        condition_chemistry = (_dirac_switch_react[_qp] > 1. ? true : false);
        cutoff = _time_react[_qp];
    }else{
        condition_chemistry = (_dirac_switch_react[_qp] >= _switch_react ? true : false);
        cutoff = 1.;
    }

    //set to zero 
    _r1[_qp] = 0.;
    _r2[_qp] = 0.;

    //add only after takeover is done
    if(_dirac_switch_react[_qp] > cutoff){ //turn on after reaction heat has gone by
        //compute the rates of reaction of all reactions
        _r1[_qp] = _Z1 * std::exp(- _E1 / (R_RDX * _T[_qp]));
        _r2[_qp] = _Z2 * std::exp(- _E2 / (R_RDX * _T[_qp]));
    }

    _Q1[_qp] = (_a1) + (_b1 * std::max(0., _T[_qp] - _T_trans));
    _Q2[_qp] = (_a2) + (_b2 * std::max(0., _T[_qp] - _T_trans));

    //compute heat from decomposition

    _Y1_dot[_qp] = - _r1[_qp] * _Y1[_qp]; //rate of change Y1
    _Y2_dot[_qp] = _r1[_qp] * _Y1[_qp] - _r2[_qp] * _Y2[_qp]; //rate of change Y2
    _Y3_dot[_qp] = _r2[_qp] * _Y2[_qp]; //rate of change Y3

    if(_use_lump){ //divide by rho*Cp so that MassLumped has consistent units
        _Y1_dot[_qp] *= 1. / (_rho[_qp] * _cv[_qp]);
        _Y2_dot[_qp] *= 1. / (_rho[_qp] * _cv[_qp]);
        _Y3_dot[_qp] *= 1. / (_rho[_qp] * _cv[_qp]);
    }
    
    _q_decomposition[_qp] = _rho[_qp] * (- _Q1[_qp] * _Y1_dot[_qp] + _Q2[_qp] * _Y3_dot[_qp]);
}