
#include "ComputeYdots.h"
#include "RankTwoTensor.h"

registerMooseObject("catApp", ComputeYdots);

InputParameters
ComputeYdots::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute the rate of reaction of each species");
    params.addRequiredCoupledVar("temperature", "temperature");
    params.addRequiredCoupledVar("Y1", "mass fraction 1");
    params.addRequiredCoupledVar("Y2", "mass fraction 2");
    params.addRequiredCoupledVar("Y3", "mass fraction 3");
    params.addRequiredCoupledVar("Y4", "mass fraction 4");
    params.addRequiredParam<Real>("Z1", "pre-exponential factor Z1");
    params.addRequiredParam<Real>("Z2", "pre-exponential factor Z2");
    params.addRequiredParam<Real>("Z3", "pre-exponential factor Z3");
    params.addRequiredParam<Real>("Rg", "thermodynamic constan");
    params.addRequiredParam<Real>("E1", "reaction energy 1");
    params.addRequiredParam<Real>("E2", "reaction energy 2");
    params.addRequiredParam<Real>("E3", "reaction energy 3");
    return params;
}

ComputeYdots::ComputeYdots(const InputParameters & parameters)
  : Material(parameters),
    _T(coupledValue("temperature")),
    _Y1(coupledValue("Y1")),
    _Y2(coupledValue("Y2")),
    _Y3(coupledValue("Y3")),
    _Y4(coupledValue("Y4")),
    //reaction rate parameters
    _Z1(getParam<Real>("Z1")),
    _Z2(getParam<Real>("Z2")),
    _Z3(getParam<Real>("Z3")),
    
    _E1(getParam<Real>("E1")),
    _E2(getParam<Real>("E2")),
    _E3(getParam<Real>("E3")),

    //declare the material properties that will store the rates
    _Y1dot(declareProperty<Real>("Y1dot")),
    _Y2dot(declareProperty<Real>("Y2dot")),
    _Y3dot(declareProperty<Real>("Y3dot")),
    _Y4dot(declareProperty<Real>("Y4dot")),
    _Rg(getParam<Real>("Rg")),

    //declare rate of reaction terms
    _r1(declareProperty<Real>("r1")),
    _dr1dT(declareProperty<Real>("dr1dT")),
    _r2(declareProperty<Real>("r2")),
    _dr2dT(declareProperty<Real>("dr2dT")),
    _r3(declareProperty<Real>("r3")),
    _dr3dT(declareProperty<Real>("dr3dT"))
{   
}

void
ComputeYdots::computeQpProperties()
{
    //compute the rates of reaction of all reactions
    _r1[_qp] = _Z1 * std::exp(- _E1 / (_Rg * _T[_qp]));
    _dr1dT[_qp] = _r1[_qp] * (_E1) * (1.0 / (_Rg * std::pow(_T[_qp], 2.0)));
    _r2[_qp] = _Z2 * std::exp(- _E2 / (_Rg * _T[_qp]));
    _dr2dT[_qp] = _r2[_qp] * (_E2) * (1.0 / (_Rg * std::pow(_T[_qp], 2.0)));
    _r3[_qp] = _Z3 * std::exp(- _E3 / (_Rg * _T[_qp]));
    _dr3dT[_qp] = _r3[_qp] * (_E3) * (1.0 / (_Rg * std::pow(_T[_qp], 2.0)));

    //use tarver 3-step reaction model to compute the Y_dot term on the 
    //RHS of the conservation of chemical species equation
    //this computes the RHS as is, with the respective sign
    _Y1dot[_qp] = - 1.0 * _r1[_qp] * _Y1[_qp];
    _Y2dot[_qp] = (_r1[_qp] * _Y1[_qp]) - (_r2[_qp] * _Y2[_qp]);
    _Y3dot[_qp] = (_r2[_qp] * _Y2[_qp]) - (_r3[_qp] * std::pow(_Y3[_qp], 2.0));
    _Y4dot[_qp] = _r3[_qp] * std::pow(_Y3[_qp], 2.0);
}
