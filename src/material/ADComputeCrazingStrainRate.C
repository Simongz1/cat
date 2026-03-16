#include "ADComputeCrazingStrainRate.h"

registerMooseObject("mlApp", ADComputeCrazingStrainRate);

InputParameters
ADComputeCrazingStrainRate::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADMaterial>::validParams();
  //params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("crazing strain rate constitutive model");

  //variables

  //parameters
  
  return params;
}

ADComputeCrazingStrainRate::ADComputeCrazingStrainRate(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ADMaterial>(parameters),
    //ADSingleVariableReturnMappingSolution(parameters),
    _gamma_c_dot(getADMaterialProperty<Real>("gamma_c_dot")),
    _switch(getADMaterialProperty<Real>("switch")),
    _nMnM(getADMaterialProperty<RankTwoTensor>("nMnM")),

    //declarations
    _epsilon_c_dot(declareADProperty<RankTwoTensor>("epsilon_c_dot"))
{}

void
ADComputeCrazingStrainRate::initialSetup()
{}

void
ADComputeCrazingStrainRate::initQpStatefulProperties()
{
  ADMaterial::initQpStatefulProperties();
  _epsilon_c_dot[_qp].zero();
}

void
ADComputeCrazingStrainRate::computeQpProperties()
{
  //condition for flow
  const Real on = (_switch[_qp] >= 0. ? 1. : 0.);
  _epsilon_c_dot[_qp] = on * _gamma_c_dot[_qp] * _nMnM[_qp];
}