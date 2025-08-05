#include "ADY4_dot.h"

registerMooseObject("beaverApp", ADY4_dot);

InputParameters
ADY4_dot::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("computes the RHS part of the conservation equation using Tarver model");
  return params;
}

ADY4_dot::ADY4_dot(const InputParameters & parameters)
  : ADKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _Y4dot(getADMaterialProperty<Real>("Y4dot"))
{}

ADReal
ADY4_dot::computeQpResidual()
{
  //computes the RHS Ydot term for the chemical species conservation
  //equaton
  //different from dY_dt computation, this directly references a material prop

  ADReal res = 0.0;
  res += - _Y4dot[_qp] * _test[_i][_qp]; //Y1_dot
  return res;
}
