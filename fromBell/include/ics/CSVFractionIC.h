#pragma once

#include "InitialCondition.h"

class CSVFractionIC : public InitialCondition
{
public:
  static InputParameters validParams();
  CSVFractionIC(const InputParameters & parameters);

protected:
  virtual Real value(const Point & p) override;

  const VariableName _fraction_var; // Aux variable name
};
