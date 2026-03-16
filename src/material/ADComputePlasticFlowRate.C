#include "ADComputePlasticFlowRate.h"

registerMooseObject("mlApp", ADComputePlasticFlowRate);

InputParameters
ADComputePlasticFlowRate::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADMaterial>::validParams();
  //params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("Plastic flow rate constitutive model");

  //variables
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredCoupledVar("shear_strength", "shear_strength");

  //parameters
  params.addRequiredParam<Real>("A", "A");
  params.addRequiredParam<Real>("gamma_p_0", "gamma_p_0");
  params.addRequiredParam<Real>("alpha_shear", "alpha_shear");

  params.addRequiredParam<Real>("h", "h");
  params.addRequiredParam<Real>("ss", "ss");
  
  return params;
}

ADComputePlasticFlowRate::ADComputePlasticFlowRate(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ADMaterial>(parameters),
    //ADSingleVariableReturnMappingSolution(parameters),
    _temperature(adCoupledValue("temperature")),
    _shear_strength(adCoupledValue("shear_strength")),
    _gamma_p_0(getParam<Real>("gamma_p_0")),
    _alpha_shear(getParam<Real>("alpha_shear")),
    _h(getParam<Real>("h")),
    _ss(getParam<Real>("ss")),

    _A(getParam<Real>("A")),
    _sp(getADMaterialProperty<RankTwoTensor>("sp")),
    _s_pressure(getADMaterialProperty<Real>("s_pressure")),

    //declarations
    _gamma_p_dot(declareADProperty<Real>("gamma_p_dot")),
    _shear_strength_dot(declareADProperty<Real>("shear_strength_dot"))
{}

void
ADComputePlasticFlowRate::initialSetup()
{
  _gamma_p_dot[_qp] = 0.;
  _shear_strength_dot[_qp] = 0.;
}

void
ADComputePlasticFlowRate::initQpStatefulProperties()
{
  ADMaterial::initQpStatefulProperties();
  _gamma_p_dot[_qp] = 0.;
  _shear_strength_dot[_qp] = 0.;
}

void
ADComputePlasticFlowRate::computeQpProperties()
{
  //compute tau
  ADRankTwoTensor sp_dev = _sp[_qp].deviatoric();
  ADReal inner_tau = sp_dev.doubleContraction(sp_dev / 2.);
  ADReal tau = MetaPhysicL::sqrt(inner_tau);

  //compute s_tilde
  ADReal s_tilde = _shear_strength[_qp] + _alpha_shear * _s_pressure[_qp];

  ADReal inner_expression;
  inner_expression = - _A * _shear_strength[_qp] / _temperature[_qp];
  inner_expression *= (1. - MetaPhysicL::pow( tau / s_tilde , 5. / 6.));
  _gamma_p_dot[_qp] = _gamma_p_0 * MetaPhysicL::exp(inner_expression);

  //compute s_dot 
  _shear_strength_dot[_qp] = _h * (1. - (_shear_strength[_qp] / _ss)) * _gamma_p_dot[_qp];
}