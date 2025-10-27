#include "CSVFractionIC.h"

registerMooseObject("mistApp", CSVFractionIC);

InputParameters
CSVFractionIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addClassDescription("read fraction from the variable storing the csv data");
  params.addRequiredCoupledVar("fraction_csv", "fraction_csv");
  return params;
}

CSVFractionIC::CSVFractionIC(const InputParameters & parameters)
  : InitialCondition(parameters),
  _fraction_csv(coupledValue("fraction_csv"))
{
}

Real
CSVFractionIC::value(const Point & p)
{
  return _fraction_csv[_qp];
}