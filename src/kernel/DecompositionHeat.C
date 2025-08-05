#include "DecompositionHeat.h"

registerMooseObject("beaverApp", DecompositionHeat);

InputParameters
DecompositionHeat::validParams()
{
  InputParameters params = HeatSource::validParams(); //we use AD to avoid complicated Jacobian computation
  params.addClassDescription("compute the heat of reaction contribution to the heat source");
  params.addRequiredParam<std::vector<Real>>("Qs", "heat of reaction values");
  return params;
}

DecompositionHeat::DecompositionHeat(const InputParameters & parameters)
  : HeatSource(parameters),
    _Qs(getParam<std::vector<Real>>("Qs")),
    _Ydot(getMaterialProperty<std::vector<Real>>("Ydot"))
{}

Real
DecompositionHeat::computeQpResidual()
{
  //NOTE: the residual equation does not directly depend on temperature
  //but the property values called from ComputeYRates do depend on temperature
  //so that the coupling is both ways.
  //define the source term value
  Real q_dec_tot = 0.0;
  int ProblemSize = _Qs.size();
  for (unsigned int i; i < ProblemSize; i++)
  {
    q_dec_tot += _Qs[i] * _Ydot[_qp][i];
  }
  return q_dec_tot * _test[_i][_qp];
}