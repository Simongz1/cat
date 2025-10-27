#pragma once

#include "AuxKernel.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "MooseUtils.h"
#include "MooseMesh.h"
#include "FEProblemBase.h"
#include "libmesh/mesh_base.h"
#include "libmesh/bounding_box.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/point.h"
#include "MooseRandom.h"
#include "SubProblem.h"
#include "SystemBase.h"

/**
 *
 */
class PolycrystalFractionAux : public AuxKernel
{
public:
  static InputParameters validParams();

  PolycrystalFractionAux(const InputParameters & parameters);
  const std::string _csv_fraction;
  const VariableValue &_density_i;
  std::vector<std::vector<Real>> _csv_total_fraction;
protected:
  virtual Real computeValue() override;
  virtual std::vector<std::vector<Real>> readCSV(const std::string csv_file_name);
};