#include "MooseUtils.h"
#include "MooseMesh.h"
#include "FEProblemBase.h"
#include "libmesh/mesh_base.h"
#include "libmesh/bounding_box.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/point.h"
#include "GeneralUserObject.h"
#include "MooseRandom.h"
#include "AuxiliarySystem.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"
#include "MooseApp.h"
#include "libmesh/dof_map.h"

class PolycrystalDensityUO : public GeneralUserObject
{
public:
  static InputParameters validParams();
  PolycrystalDensityUO(const InputParameters & parameters);
  virtual void execute() override;  
  virtual void initialize() override;
  virtual void finalize() override; 
  virtual void initialSetup() override;
  virtual std::vector<std::vector<Real>> readCSV(const std::string csv_file_name);
protected:
  unsigned int _num_grains;
  const std::vector<unsigned int> _target_grains;
  const std::vector<unsigned int> _range_in;
  const std::vector<unsigned int> _range_out;
  const bool _generate_matrix;
  const Real _max_grain_size;
  const Real _min_center_spacing;
  const Real _matrix_thickness;
  const std::vector<unsigned int> _sizes;
  const std::vector<Real> _sizes_fraction;
  const std::string _csv_fraction;
  //const VariableValue &_density_i;

  ////////////////
  std::vector<std::vector<Real>> _csv_total_fractions;
  
  std::vector<Point> _centers;
  std::vector<Real> _radii;

};
