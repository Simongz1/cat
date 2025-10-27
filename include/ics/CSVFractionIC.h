#pragma once

#include "InitialCondition.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

/**
 *
 */
class CSVFractionIC : public InitialCondition
{
public:
  static InputParameters validParams();

  CSVFractionIC(const InputParameters & parameters);
  const VariableValue &_fraction_csv;

  /////////////////////////////////
  std::vector<std::vector<Real>> _csv_total_density;
  virtual Real value(const Point & p) override;
};