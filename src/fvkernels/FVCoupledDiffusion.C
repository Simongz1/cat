#include "FVCoupledDiffusion.h"
#include "FVDiffusion.h"

registerADMooseObject("mlApp", FVCoupledDiffusion);

InputParameters
FVCoupledDiffusion::validParams()
{
  InputParameters params = FVDiffusion::validParams();
  params.addClassDescription("returns the advection of a coupled variable as residual contribution");
  params.addRequiredCoupledVar("phi", "name of the coupled variable");
  params.addParam<Real>("factor", 1.0, "factor that premultiplies the residual");
  return params;
}

FVCoupledDiffusion::FVCoupledDiffusion(const InputParameters & parameters)
  : FVDiffusion(parameters),
    _coupled_variable(dynamic_cast<const MooseVariableFV<Real> &>(*getFieldVar("phi", 0))),
    _factor(getParam<Real>("factor"))
{}

ADReal
FVCoupledDiffusion::computeQpResidual()
{ 
  //obtain the required gradient times normal
  ADReal dvardnormal = Moose::FV::gradUDotNormal(*_face_info, _coupled_variable, determineState(), _correct_skewness);

  //return statement
  return _factor * dvardnormal;
}