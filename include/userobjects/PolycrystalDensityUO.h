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
  virtual Real assignPoreValue(unsigned int grain_id, const Elem * elem, bool is_grain, Real rand_value);
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
  const std::string _csv_fraction_pore;
  const bool _bulk_grains;
  const unsigned int _n_cracks;
  const unsigned int _l_cracks;
  const unsigned int _n_pores;

  const unsigned int _bulk_MicroID;
  const Real _bulk_RDX_fraction;
  const std::vector<unsigned int> _range_pore;
  const Real _pore_RDX_fraction;


  ////////////////
  std::vector<std::vector<Real>> _csv_total_fractions;
  std::vector<std::vector<Real>> _csv_total_fractions_pore;
  
  std::vector<Point> _centers;
  std::vector<Real> _radii;
  std::vector<Real> _grainID;

};
