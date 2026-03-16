#include "ADComputePlasticStrainRate.h"

registerMooseObject("mlApp", ADComputePlasticStrainRate);

InputParameters
ADComputePlasticStrainRate::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADMaterial>::validParams();
  //params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("plastic strain rate constitutive model");

  //variables

  //parameters
  
  return params;
}

ADComputePlasticStrainRate::ADComputePlasticStrainRate(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ADMaterial>(parameters),
    //ADSingleVariableReturnMappingSolution(parameters),
    _gamma_p_dot(getADMaterialProperty<Real>("gamma_p_dot")),
    _sp(getADMaterialProperty<RankTwoTensor>("sp")),
    _switch(getADMaterialProperty<Real>("switch")),

    //declarations
    _epsilon_p_dot(declareADProperty<RankTwoTensor>("epsilon_p_dot"))
{}

void
ADComputePlasticStrainRate::initialSetup()
{}

void
ADComputePlasticStrainRate::initQpStatefulProperties()
{
  ADMaterial::initQpStatefulProperties();
  _epsilon_p_dot[_qp].zero();
}

void
ADComputePlasticStrainRate::computeQpProperties()
{
  //condition for flow
  const Real on = (_switch[_qp] >= 0. ? 1. : 0.);
  
  //compute deviatoric small plastic stress
  ADRankTwoTensor sp_dev = _sp[_qp].deviatoric();
  ADReal sp_dev_norm = MetaPhysicL::sqrt(sp_dev.doubleContraction(sp_dev));
  _epsilon_p_dot[_qp] = on * _gamma_p_dot[_qp] * sp_dev / sp_dev_norm;
}