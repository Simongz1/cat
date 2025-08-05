#include "IGdotRDX.h"

registerMooseObject("catApp", IGdotRDX);

InputParameters
IGdotRDX::validParams()
{
  InputParameters params = Kernel::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute the heat of reaction contribution to the heat source");
  params.addRequiredParam<Real>("x0", "x0");
  params.addRequiredParam<Real>("y0", "x0");
  params.addRequiredParam<Real>("a", "a");

  params.addRequiredParam<Real>("x1", "x1");
  params.addRequiredParam<Real>("y1", "y1");
  params.addRequiredParam<Real>("z1", "z1");

  params.addRequiredParam<Real>("x2", "x2");
  params.addRequiredParam<Real>("y2", "y2");
  params.addRequiredParam<Real>("z2", "z2");

  params.addRequiredParam<Real>("I", "I");
  params.addRequiredParam<Real>("G1", "G1");
  params.addRequiredParam<Real>("G2", "G2");  

  params.addRequiredParam<Real>("limI", "limI");
  params.addRequiredParam<Real>("limG1", "limG1");
  params.addRequiredCoupledVar("lambda", "lambda");
  return params;
}

IGdotRDX::IGdotRDX(const InputParameters & parameters)
  : Kernel(parameters),
    _x0(getParam<Real>("x0")),
    _y0(getParam<Real>("y0")),
    _a(getParam<Real>("a")),

    _x1(getParam<Real>("x1")),
    _y1(getParam<Real>("y1")),
    _z1(getParam<Real>("z1")),

    _x2(getParam<Real>("x2")),
    _y2(getParam<Real>("y2")),
    _z2(getParam<Real>("z2")),

    _I(getParam<Real>("I")),
    _G1(getParam<Real>("G1")),
    _G2(getParam<Real>("G2")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _pressure_total(getMaterialProperty<Real>("pressure_total")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),

    _limI(getParam<Real>("limI")),
    _limG1(getParam<Real>("limG1")),
    _lambda(coupledValue("lambda"))
{}

Real
IGdotRDX::computeQpResidual()
{
  //define rate terms
  RankTwoTensor F = _deformation_gradient[_qp];
  Real J = F.det();
  Real P = std::abs(_pressure_total[_qp]);
  Real deltaV = 1. / (1. + _mechanical_strain[_qp].trace());

  Real switch1 = (_u[_qp] < _limI) ? 1. : 0.;
  Real switch2 = (_u[_qp] < _limG1) ? 1. : 0.;
  Real switch3 = (_u[_qp] > _limG1 && _u[_qp] < 1.) ? 1. : 0.;

  Real IG = _I * std::pow((1. - _lambda[_qp]), _x0) * std::pow(((deltaV) - 1 - _a), _y0) * switch1;
  Real G1 = _G1 * std::pow((1. - _lambda[_qp]), _x1) * std::pow(_lambda[_qp], _y1) * std::pow(P, _z1) * switch2;
  Real G2 = _G2 * std::pow((1. - _lambda[_qp]), _x2) * std::pow(_lambda[_qp], _y2) * std::pow(P, _z2) * switch3;

  Real rate = IG + G1 + G2;
  return - rate * _test[_i][_qp];
}

Real
IGdotRDX::computeQpJacobian()
{
  RankTwoTensor F = _deformation_gradient[_qp];
  Real J = F.det();
  Real P = std::abs(_pressure_total[_qp]);
  Real deltaV = 1. / (1. + _mechanical_strain[_qp].trace());

  Real switch1 = (_u[_qp] < _limI) ? 1. : 0.;
  Real dIG = _I * (- _x0) * std::pow((1. - _u[_qp]), _x0 - 1.) * std::pow(((deltaV) - 1 - _a), _y0) * switch1;

  Real switch2 = (_u[_qp] < _limG1) ? 1. : 0.;
  Real dG1 = - _G1 * _x1 * std::pow((1. - _u[_qp]), _x1 - 1.) * std::pow(_u[_qp], _y1) * std::pow(P, _z1)
            + _G1 * std::pow((1. - _u[_qp]), _x1) * _y1 * std::pow(_u[_qp], _y1 - 1.) * std::pow(P, _z1);
  dG1 *= switch2;

  Real switch3 = (_u[_qp] > _limG1 && _u[_qp] < 1.) ? 1. : 0.;
  Real dG2 = - _G2 * _x2 * std::pow((1. - _u[_qp]), _x2 - 1.) * std::pow(_u[_qp], _y2) * std::pow(P, _z2)
            + _G2 * std::pow((1. - _u[_qp]), _x2) * _y2 * std::pow(_u[_qp], _y2 - 1.) * std::pow(P, _z2);
  dG2 *= switch3;

  Real Jac = dIG + dG1 + dG2;
  return - Jac * _test[_i][_qp] * _phi[_j][_qp];
}

// Real
// IGdotRDX::computeQpOffDiagonalJacobian(unsigned int jvar)
// {
//   return - 0.0 * _test[_i][_qp] * _phi[_j][_qp];
// }