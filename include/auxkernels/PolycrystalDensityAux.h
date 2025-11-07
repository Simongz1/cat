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
  const std::string _csv_density_pore;
  const VariableValue &_density_i;
  const Real _density_scaling;

  //parameters for bulk grains
  const unsigned int _bulk_MicroID;
  const Real _bulk_RDX_density;
  const std::vector<unsigned int> _range_pore;
  std::vector<std::vector<Real>> _csv_total_density;
  std::vector<std::vector<Real>> _csv_total_density_pore;
protected:
  virtual Real computeValue() override;
  virtual std::vector<std::vector<Real>> readCSV(const std::string csv_file_name);
};