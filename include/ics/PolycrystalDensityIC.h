#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "RandomIC.h"
#include "MooseUtils.h"
#include "MooseMesh.h"
#include "FEProblemBase.h"
#include "libmesh/mesh_base.h"
#include "libmesh/bounding_box.h"
#include "libmesh/mesh_tools.h"

/**
 *
 */
class PolycrystalDensityIC : public RandomIC
{
public:
  static InputParameters validParams();

  PolycrystalDensityIC(const InputParameters & parameters);
  const Real _num_grains;
  const Real _min_in;
  const Real _max_in;
  const Real _min_out;
  const Real _max_out;
  const Real _target_grain;
  std::vector<Point> _centers;

  /////////////////////////////////
  virtual Real value(const Point & p) override;
  virtual void initialSetup() override;
};