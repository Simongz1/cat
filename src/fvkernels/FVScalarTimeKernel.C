#include "FVScalarTimeKernel.h"
#include "NS.h"
#include "SystemBase.h"

registerADMooseObject("mlApp", FVScalarTimeKernel);

InputParameters
FVScalarTimeKernel::validParams()
{
  InputParameters params = FVTimeKernel::validParams();
  params.addClassDescription("Residual contribution from time derivative of a variable for the finite volume method times a scalar functor.");
  params.addRequiredParam<MooseFunctorName>("scalar", "name of the scalar quantity");
  return params;
}

FVScalarTimeKernel::FVScalarTimeKernel(const InputParameters & parameters)
  : FVTimeKernel(parameters),
    _scalar(getFunctor<ADReal>("scalar"))
{
}

ADReal
FVScalarTimeKernel::computeQpResidual()
{
  //form element argument
  const auto elem_arg = makeElemArg(_current_elem);
  const ADReal scalar = _scalar(elem_arg, determineState());
  return scalar * _u_dot[_qp];
}