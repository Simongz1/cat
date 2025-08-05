#include "ADY2_dot_RDX.h"

registerMooseObject("catApp", ADY2_dot_RDX);

InputParameters
ADY2_dot_RDX::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Y1 species time derivate term of the conservation equation"
                              "this kernel computes only the time derivative contribution"
                              "to the residual");
  params.addRequiredCoupledVar("Y1", "Y1");
  return params;
}

ADY2_dot_RDX::ADY2_dot_RDX(const InputParameters & parameters)
  : ADKernel(parameters),
    _r2(getADMaterialProperty<Real>("r2")),
    _r1(getADMaterialProperty<Real>("r1")),
    _Y1(coupledValue("Y1"))
{}

ADReal
ADY2_dot_RDX::computeQpResidual()
{
  ADReal res;
  res = (_r1[_qp] * _Y1[_qp]) - (_r2[_qp] * _u[_qp]); //Y2_dot
  return - res * _test[_i][_qp];
}