#include "ADSurrogateY.h"

registerMooseObject("catApp", ADSurrogateY);

InputParameters
ADSurrogateY::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("append the rate of surrogate evolution into the residual of each chemical variable");
  params.addRequiredParam<MaterialPropertyName>("surrogate_rate_name", "the name of the surrogate rate");
  params.addRequiredParam<bool>("positive", "whether the residual is positive or negative");
  return params;
}

ADSurrogateY::ADSurrogateY(const InputParameters & parameters)
  : ADKernel(parameters),
    _surrogate_name(getParam<MaterialPropertyName>("surrogate_rate_name")),
    _surrogate_rate(getADMaterialPropertyByName<Real>(_surrogate_name)),
    _positive(getParam<bool>("positive"))
{}

ADReal
ADSurrogateY::computeQpResidual()
{
  Real sign;
  if (_positive){
    sign = 1.;
  }else{
    sign = - 1.;
  }
  return sign * _surrogate_rate[_qp] * _test[_i][_qp];
}