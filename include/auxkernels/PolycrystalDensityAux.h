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
#include "PolycrystalDensityUO.h"
#include "SubProblem.h"
#include "SystemBase.h"

/**
 *
 */
class PolycrystalDensityAux : public AuxKernel
{
public:
  static InputParameters validParams();

  PolycrystalDensityAux(const InputParameters & parameters);
  const std::string _csv_density;
  const VariableValue &_density_i;
  const Real _density_scaling;
  std::vector<std::vector<Real>> _csv_total_density;
protected:
  virtual Real computeValue() override;
  virtual std::vector<std::vector<Real>> readCSV(const std::string csv_file_name);
};