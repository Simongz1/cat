#include "ADComputePolymerStress.h"

registerMooseObject("mlApp", ADComputePolymerStress);

InputParameters
ADComputePolymerStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes the compressible Neo-Hookean stress for the polymer phase.");
  params.addParam<MaterialPropertyName>("bulk_modulus", "bulk_modulus", "name of the bulk modulus property");
  params.addParam<MaterialPropertyName>("poisson_ratio", "poisson_ratio", "name of the poisson ratio property");
  return params;
}

ADComputePolymerStress::ADComputePolymerStress(const InputParameters & parameters)
  : Material(parameters),
    //obtain bulk modulus from name
    _bulk_modulus_name(getParam<MaterialPropertyName>("bulk_modulus")),
    _bulk_modulus(getADMaterialProperty<Real>(_bulk_modulus_name)),

    //obtain poisson ratio from name
    _poisson_ratio_name(getParam<MaterialPropertyName>("poisson_ratio")),
    _poisson_ratio(getADMaterialProperty<Real>(_poisson_ratio_name)),

    //obtain kinematic variables
    _be_bar(getADMaterialProperty<RankTwoTensor>("be_bar")),
    _J(getADMaterialProperty<Real>("J")),

    //declare the binder stress
    _sigma_binder(declareADProperty<RankTwoTensor>("sigma_binder"))
{}

void
ADComputePolymerStress::initialSetup()
{}

void
ADComputePolymerStress::initQpStatefulProperties()
{}

void
ADComputePolymerStress::computeQpProperties()
{ 
  ADRankTwoTensor I2;
  I2.setToIdentity();

  //compute shear modulus from bulk and poisson ratio
  ADReal shear_modulus = 3.0 * _bulk_modulus[_qp] * (1.0 - 2.0 * _poisson_ratio[_qp]);
  shear_modulus *= 1.0 / (2.0 * (1.0 + _poisson_ratio[_qp]));

  //compute stress
  _sigma_binder[_qp] = _bulk_modulus[_qp] * (_J[_qp] - 1.0) * I2;
  _sigma_binder[_qp] += (1.0 / _J[_qp]) * shear_modulus * _be_bar[_qp].deviatoric();
}