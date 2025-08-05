#include "ComputeYdotsRDXNS.h"
#include "RankTwoTensor.h"

registerMooseObject("catApp", ComputeYdotsRDXNS);

InputParameters
ComputeYdotsRDXNS::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute the rate of reaction of each species for a 3 phase - 2 stage reaction model");
    params.addRequiredCoupledVar("temperature", "temperature");
    params.addRequiredCoupledVar("Y1", "mass fraction 1");
    params.addRequiredCoupledVar("Y2", "mass fraction 2");
    params.addRequiredCoupledVar("Y3", "mass fraction 3");
    params.addRequiredParam<Real>("Z1", "pre-exponential factor Z1");
    params.addRequiredParam<Real>("Z2", "pre-exponential factor Z2");
    params.addRequiredParam<Real>("Rg", "thermodynamic constan");
    params.addRequiredParam<Real>("MW", "Molecular Weight of RDX");
    params.addRequiredParam<Real>("E1", "reaction energy 1");
    params.addRequiredParam<Real>("E2", "reaction energy 2");
    params.addRequiredParam<Real>("T_trans", "transition temperature for heat of reaction");
    params.addRequiredParam<Real>("a1", "a1 parameter for heat of reaction 1");
    params.addRequiredParam<Real>("b1", "b1 parameter for heat of reaction 1");
    params.addRequiredParam<Real>("a2", "a2 parameter for heat of reaction 2");
    params.addRequiredParam<Real>("b2", "b2 parameter for heat of reaction 2");
    params.addCoupledVar("vx", "vx");
    params.addRequiredCoupledVar("dirac_switch_react", "dirac");
    params.addRequiredParam<Real>("rate_limit","reaction rate limit");
    params.addRequiredParam<Real>("switch_react", "dirac switch value to turn on reaction");

    return params;
}

ComputeYdotsRDXNS::ComputeYdotsRDXNS(const InputParameters & parameters)
  : Material(parameters),
    _T(coupledValue("temperature")),
    _Y1(coupledValue("Y1")),
    _Y2(coupledValue("Y2")),
    _Y3(coupledValue("Y3")),
    //reaction rate parameters
    _Z1(getParam<Real>("Z1")),
    _Z2(getParam<Real>("Z2")),
    
    _E1(getParam<Real>("E1")),
    _E2(getParam<Real>("E2")),

    //declare the material properties that will store the rates
    _Y1dot(declareProperty<Real>("Y1dot")),
    _Y2dot(declareProperty<Real>("Y2dot")),
    _Y3dot(declareProperty<Real>("Y3dot")),
    _Rg(getParam<Real>("Rg")),
    _MW(getParam<Real>("MW")),

    //declare rate of reaction terms
    _r1(declareProperty<Real>("r1")),
    _dr1dT(declareProperty<Real>("dr1dT")),
    _r2(declareProperty<Real>("r2")),
    _dr2dT(declareProperty<Real>("dr2dT")),

    //declare heats of reaction
    _Q1(declareProperty<Real>("Q1")),
    _Q2(declareProperty<Real>("Q2")),
    _a1(getParam<Real>("a1")),
    _b1(getParam<Real>("b1")),
    _a2(getParam<Real>("a2")),
    _b2(getParam<Real>("b2")),
    _T_trans(getParam<Real>("T_trans")),
    _dirac_switch_react(coupledValue("dirac_switch_react")),
    //rate limit
    _rate_limit(getParam<Real>("rate_limit")),
    _switch_react(getParam<Real>("switch_react"))
{   
}

void
ComputeYdotsRDXNS::computeQpProperties()
{
    Real R_RDX = _Rg;
    if(_dirac_switch_react[_qp] >= _switch_react){ //turn on after reaction heat has gone by
        //compute the rates of reaction of all reactions
        _r1[_qp] = _Z1 * std::exp(- _E1 / (R_RDX * _T[_qp]));
        _r1[_qp] = std::min(_r1[_qp], _rate_limit); //limit the rate of reaction
        _dr1dT[_qp] = _r1[_qp] * (_E1) * (1.0 / (R_RDX * std::pow(_T[_qp], 2.0)));
        _r2[_qp] = _Z2 * std::exp(- _E2 / (R_RDX * _T[_qp]));
        _r2[_qp] = std::min(_r2[_qp], _rate_limit); //limit the rate of reaction
        _dr2dT[_qp] = _r2[_qp] * (_E2) * (1.0 / (R_RDX * std::pow(_T[_qp], 2.0)));
    }else{ //don't do anything before
        _r1[_qp] = 0.0;
        _dr1dT[_qp] = _r1[_qp] * (_E1) * (1.0 / (R_RDX * std::pow(_T[_qp], 2.0)));
        _r2[_qp] = 0.0;
        _dr2dT[_qp] = _r2[_qp] * (_E2) * (1.0 / (R_RDX * std::pow(_T[_qp], 2.0)));
    }


    //use tarver 3-step reaction model to compute the Y_dot term on the 
    //RHS of the conservation of chemical species equation
    //this computes the RHS as is, with the respective sign
    _Y1dot[_qp] = - 1.0 * _r1[_qp] * _Y1[_qp];
    _Y2dot[_qp] = (_r1[_qp] * _Y1[_qp]) - (_r2[_qp] * _Y2[_qp]);
    _Y3dot[_qp] = (_r2[_qp] * _Y2[_qp]);


    //compute the heat of reaction for RDX
    Real Q1;
    Real Q2;
    // multiplies with r1 r2 in kernels
    _Q1[_qp] = (_a1) + _b1 * std::max(0., _T[_qp] - _T_trans);
    _Q2[_qp] = (_a2) + _b2 * std::max(0., _T[_qp] - _T_trans);  
}