#include "ADY2_dot.h"

registerMooseObject("beaverApp", ADY2_dot);

InputParameters
ADY2_dot::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("computes the RHS part of the conservation equation using Tarver model");
  return params;
}

ADY2_dot::ADY2_dot(const InputParameters & parameters)
  : ADKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _Y2dot(getADMaterialProperty<Real>("Y2dot"))
{}

ADReal
ADY2_dot::computeQpResidual()
{
  //computes the RHS Ydot term for the chemical species conservation
  //equaton
  //different from dY_dt computation, this directly references a material prop

  ADReal res = 0.0;
  res += - _Y2dot[_qp] * _test[_i][_qp]; //Y1_dot
  return res;
}
