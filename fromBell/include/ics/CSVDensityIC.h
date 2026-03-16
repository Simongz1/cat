#pragma once

#include "InitialCondition.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

/**
 *
 */
class CSVDensityIC : public InitialCondition
{
public:
  static InputParameters validParams();

  CSVDensityIC(const InputParameters & parameters);
  const std::string _csv_density;
  const VariableValue &_density_i;
  const Real _density_scaling;

  /////////////////////////////////
  std::vector<std::vector<Real>> _csv_total_density;
  virtual Real value(const Point & p) override;
  virtual std::vector<std::vector<Real>> readCSV(const std::string csv_name);
};