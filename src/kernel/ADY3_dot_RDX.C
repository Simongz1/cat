#include "ADY3_dot_RDX.h"

registerMooseObject("catApp", ADY3_dot_RDX);

InputParameters
ADY3_dot_RDX::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Y1 species time derivate term of the conservation equation"
                              "this kernel computes only the time derivative contribution"
                              "to the residual");
  params.addRequiredCoupledVar("Y2", "Y2");
  return params;
}

ADY3_dot_RDX::ADY3_dot_RDX(const InputParameters & parameters)
  : ADKernel(parameters),
    _r2(getADMaterialProperty<Real>("r2")),
    _Y2(coupledValue("Y2"))
{}

ADReal
ADY3_dot_RDX::computeQpResidual()
{
  ADReal res;
  res = (_r2[_qp] * _Y2[_qp]);
  return - res * _test[_i][_qp];
}