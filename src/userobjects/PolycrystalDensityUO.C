#include "PolycrystalDensityUO.h"

registerMooseObject("mistApp", PolycrystalDensityUO);

InputParameters
PolycrystalDensityUO::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addRequiredParam<unsigned int>("num_grains", "Number of Voronoi seeds");
  params.addRequiredParam<std::vector<unsigned int>>("target_grains", "target_grains");
  params.addRequiredParam<std::vector<unsigned int>>("range_in", "range_in");
  params.addRequiredParam<std::vector<unsigned int>>("range_out", "range_out");
  params.addRequiredParam<bool>("generate_matrix", "generate_matrix");
  params.addRequiredParam<Real>("max_grain_size", "max_grain_size");
  params.addRequiredParam<Real>("min_center_spacing", "min_center_spacing");
  params.addRequiredParam<Real>("matrix_thickness", "matrix_thickness");
  params.addRequiredParam<std::vector<unsigned int>>("sizes", "sizes");
  params.addRequiredParam<std::vector<Real>>("sizes_fraction", "sizes_fraction");

  // for fractions
  params.addRequiredParam<std::string>("csv_fraction", "csv_fraction");
  //params.addRequiredCoupledVar("density_i", "density_i");
  return params;
}

PolycrystalDensityUO::PolycrystalDensityUO(const InputParameters & params)
  : GeneralUserObject(params), 
    _num_grains(getParam<unsigned int>("num_grains")),
    _target_grains(getParam<std::vector<unsigned int>>("target_grains")),
    _range_in(getParam<std::vector<unsigned int>>("range_in")),
    _range_out(getParam<std::vector<unsigned int>>("range_out")),
    _generate_matrix(getParam<bool>("generate_matrix")),
    _max_grain_size(getParam<Real>("max_grain_size")),
    _min_center_spacing(getParam<Real>("min_center_spacing")),
    _matrix_thickness(getParam<Real>("matrix_thickness")),
    _sizes(getParam<std::vector<unsigned int>>("sizes")),
    _sizes_fraction(getParam<std::vector<Real>>("sizes_fraction")),
    _csv_fraction(getParam<std::string>("csv_fraction"))
{
    _csv_total_fractions = readCSV(_csv_fraction);
}

void
PolycrystalDensityUO::initialSetup(){
    // Generate Voronoi centers once before AuxVariables initialize
  const BoundingBox bbox = MeshTools::create_bounding_box(_fe_problem.mesh().getMesh());
  const Point min_corner = bbox.min();
  const Point max_corner = bbox.max();

  _centers.clear();
  _centers.reserve(_num_grains);

  _radii.clear();
  _radii.reserve(_num_grains);

  for (unsigned int i = 0; i < _num_grains; ++i)
  { 
    //generate minimum spacing
    bool accepted = false;
    Point candidate;

    for (unsigned int attempt = 0; attempt < 1e3 && !accepted; ++attempt){

    const Real x = std::tanh(MooseRandom::randNormal()) * (max_corner(0) - min_corner(0)) + min_corner(0);
    const Real y = std::tanh(MooseRandom::randNormal()) * (max_corner(1) - min_corner(1)) + min_corner(1);
    const Real z = (_fe_problem.mesh().dimension() == 3)
                     ? std::tanh(MooseRandom::randNormal()) * (max_corner(2) - min_corner(2)) + min_corner(2)
                     : 0.0;
    candidate = Point(x, y, z);

    accepted = true;
    for (const auto & existing : _centers){
      const Real dist = (candidate - existing).norm();
      if (dist < _min_center_spacing){
        accepted = false;
        break;
      }
    }
    if (accepted){
      _centers.push_back(candidate);

      const Real r = std::tanh(MooseRandom::randNormal());
      if (r < _sizes_fraction[0]){
        _radii.push_back(_sizes[0]);
      }else{
        _radii.push_back(_sizes[1]);
      }
    }
    }
  }

  if (_tid == 0)
    mooseInfo("Generated ", _num_grains, " Voronoi centers in PolycrystalDensityUO");

  auto & sys = _fe_problem.getAuxiliarySystem();
  auto & nl_sys = _fe_problem.getNonlinearSystem(0);
  auto & var = sys.getVariable(_tid, "density_i");
  auto & var_Y1 = nl_sys.getVariable(_tid, "Y1");
  const DofMap & dof_map = sys.system().get_dof_map();
  const DofMap & nl_dof_map = nl_sys.system().get_dof_map();

  //unordered map for density
  std::unordered_map<dof_id_type, Real> elem_density;

  // Write initial values to the variable field
  for (const auto & elem : _fe_problem.mesh().getMesh().active_element_ptr_range())
  {
    Point centroid = elem->vertex_average();

    Real min_dist = std::numeric_limits<Real>::max();
    Real second_min_dist = std::numeric_limits<Real>::max();
    unsigned int nearest = 0;

    for (unsigned int i = 0; i < _centers.size(); ++i)
    {
      const Real d = (centroid - _centers[i]).norm();
      if (d < min_dist)
      {
        second_min_dist = min_dist;
        min_dist = d;
        nearest = i;
      }
      else if (d < second_min_dist){
        second_min_dist = d;
      }
    }
    //determine if its in a grain or outside
    const Real boundary_gap = (second_min_dist - min_dist);
    bool is_boundary = boundary_gap < _matrix_thickness;
    bool is_far = min_dist > _radii[nearest];
    bool is_grain = !is_boundary && !is_far;
    //
    const bool in_target =
        std::find(_target_grains.begin(), _target_grains.end(), nearest + 1) != _target_grains.end();

    const Real rand_value = MooseRandom::rand();
    const Real density_val =
        is_grain ? (_range_in[0] + rand_value * (_range_in[1] - _range_in[0]))
                  : (_range_out[0] + rand_value * (_range_out[1] - _range_out[0]));
    
    elem_density[elem->id()] = density_val;

    std::vector<dof_id_type> dof_indices;
    dof_map.dof_indices(elem, dof_indices, var.number());

    for (auto dof : dof_indices){
      sys.solution().set(dof, density_val);
    }
  }
  sys.solution().close();

  //read csv using density_i
  const std::vector<Real> data = _csv_total_fractions[0];
  std::unordered_map<dof_id_type, std::pair<Real, unsigned int>> density_accum;

  //open nonlinear system to set fraction values
  
  for (const auto & elem : _fe_problem.mesh().getMesh().active_element_ptr_range()){
    const Real val = elem_density[elem->id()];
    for (unsigned int i = 0; i < elem->n_nodes(); ++i){
      const Node * node = elem->node_ptr(i);
      const dof_id_type node_id = node->id();
      auto & pair = density_accum[node_id];
      pair.first += val;
      pair.second += 1;
    }
  }

  std::unordered_map<dof_id_type, Real> node_density;
  node_density.reserve(density_accum.size());
  for (const auto & kv : density_accum){
    node_density[kv.first] = kv.second.first / kv.second.second;
  }

  for (const auto * node : _fe_problem.mesh().getMesh().node_ptr_range()){
    std::vector<dof_id_type> dof_indices_Y1;
    nl_dof_map.dof_indices(node, dof_indices_Y1, var_Y1.number());
    
    if (dof_indices_Y1.empty()){
      continue;
    }
    const dof_id_type nid = node->id();
    const Real density_value = node_density.count(nid) ? node_density[nid] : 0.;

    const Real predicted_fraction = data[density_value];
    nl_sys.solution().set(dof_indices_Y1[0], predicted_fraction);
  }
  nl_sys.solution().close();
}

void
PolycrystalDensityUO::execute()
{
  
}

void
PolycrystalDensityUO::initialize(){

}

void 
PolycrystalDensityUO::finalize(){

}

std::vector<std::vector<Real>>
PolycrystalDensityUO::readCSV(const std::string csv_file_name){
  std::ifstream file(csv_file_name);
  if (!file.is_open()){
    mooseError("can't open CSV file");
  }

  std::vector<std::vector<Real>> data;
  std::string line;
  while (std::getline(file, line)){
    std::vector<Real> row;
    std::stringstream ss(line);
    std::string value;

    while (std::getline(ss, value, ',')){
      try{
        row.push_back(static_cast<Real>(std::stod(value)));
      }
      catch (const std::invalid_argument &e){
        mooseError("Invalid value in CSV file: " + value);
      }
    }
    data.push_back(row);
  }
  file.close();
  return data;
}