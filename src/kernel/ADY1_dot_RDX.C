#include "ADY1_dot_RDX.h"

registerMooseObject("catApp", ADY1_dot_RDX);

InputParameters
ADY1_dot_RDX::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("compute tarver model reaction rate for Y1");
  return params;
}

ADY1_dot_RDX::ADY1_dot_RDX(const InputParameters & parameters)
  : ADKernel(parameters),
    _r1(getADMaterialProperty<Real>("r1"))
{}

ADReal
ADY1_dot_RDX::computeQpResidual()
{
  ADReal res;
  res = - _r1[_qp] * _u[_qp]; //- r1 * Y1
  return - res * _test[_i][_qp];
}