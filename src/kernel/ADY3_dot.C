#include "ADY3_dot.h"

registerMooseObject("beaverApp", ADY3_dot);

InputParameters
ADY3_dot::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("computes the RHS part of the conservation equation using Tarver model");
  return params;
}

ADY3_dot::ADY3_dot(const InputParameters & parameters)
  : ADKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _Y3dot(getADMaterialProperty<Real>("Y3dot"))
{}

ADReal
ADY3_dot::computeQpResidual()
{
  //computes the RHS Ydot term for the chemical species conservation
  //equaton
  //different from dY_dt computation, this directly references a material prop

  ADReal res = 0.0;
  res += - _Y3dot[_qp] * _test[_i][_qp]; //Y1_dot
  return res;
}
