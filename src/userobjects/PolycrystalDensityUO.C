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
  params.addRequiredParam<std::string>("csv_fraction_pore", "csv_fraction_pore");

  //additional features

  params.addRequiredParam<bool>("bulk_grains", "bulk_grains");
  params.addRequiredParam<unsigned int>("n_cracks", "n_cracks");
  params.addRequiredParam<unsigned int>("l_cracks", "l_cracks");
  params.addRequiredParam<unsigned int>("n_pores", "n_pores");

  //provide bulk RDX values
  params.addRequiredParam<unsigned int>("bulk_MicroID", "bulk_MicroID");
  params.addRequiredParam<Real>("bulk_RDX_fraction", "bulk_RDX_fraction");
  params.addRequiredParam<std::vector<unsigned int>>("range_pore", "range_pore");
  params.addRequiredParam<Real>("pore_RDX_fraction", "pore_RDX_fraction");
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
    _csv_fraction(getParam<std::string>("csv_fraction")),
    _csv_fraction_pore(getParam<std::string>("csv_fraction_pore")),
    _bulk_grains(getParam<bool>("bulk_grains")),
    _n_cracks(getParam<unsigned int>("n_cracks")),
    _l_cracks(getParam<unsigned int>("l_cracks")),
    _n_pores(getParam<unsigned int>("n_pores")),

    //bulk parameters
    _bulk_MicroID(getParam<unsigned int>("bulk_MicroID")),
    _bulk_RDX_fraction(getParam<Real>("bulk_RDX_fraction")),

    //pore parameters
    _range_pore(getParam<std::vector<unsigned int>>("range_pore")),
    _pore_RDX_fraction(getParam<Real>("pore_RDX_fraction"))
{
    _csv_total_fractions = readCSV(_csv_fraction);
    _csv_total_fractions_pore = readCSV(_csv_fraction_pore);
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
      //here we can assign a grain ID based on the nearest neighbor to each element

    }
    }
  }

  if (_tid == 0)
    mooseInfo("Generated ", _num_grains, " Voronoi centers in PolycrystalDensityUO. Next step is to assign defects inside grains and nanoPBXs at the interfaces");


  //here we perform the pore element assignment
  std::unordered_map<unsigned int, std::unordered_set<const Elem *>> pores_per_grain;

  //in case we supply number of pores

  if (_n_pores > 0){
    MooseRandom::seed(1234);
    for (unsigned int gid = 0; gid < _num_grains; ++gid){
      unsigned int count = 0;
      for (const auto & elem : _fe_problem.mesh().getMesh().active_element_ptr_range()){
        if (count >= _n_pores){
          break;
        }

        Point centroid = elem->vertex_average();
        const Real d = (centroid - _centers[gid]).norm();
        if (d < _radii[gid]){
          const Real r = MooseRandom::rand();

          if (r < 0.01){ //this threshold can be changed to generate more pore sites
            pores_per_grain[gid].insert(elem);
            count++;
          }
        }
      }
    }
  }

  //here starts the variable assignment
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
    bool is_grain = !is_boundary && !is_far; //this determines which element is inside a grain
    //
    const bool in_target =
        std::find(_target_grains.begin(), _target_grains.end(), nearest + 1) != _target_grains.end();

    const Real rand_value = MooseRandom::rand();

    //create branch for the case when _bulk_grains = true
    //this assigns an artificial/placeholder MicroID to all elements inside grains
    Real density_val;

    if (_bulk_grains){
      if (pores_per_grain.count(nearest) && pores_per_grain[nearest].count(elem)){
        const unsigned int n_types = _range_pore.size();
        const unsigned int idx = static_cast<unsigned int>(std::floor(MooseRandom::rand() * n_types)) % n_types;

        density_val = static_cast<Real>(_range_pore[idx]);
      }
      else if (is_grain){
        density_val = static_cast<Real>(_bulk_MicroID);
      }
      else{
        density_val = _range_out[0] + rand_value * (_range_out[1] - _range_out[0]);
      }
    }
    else {
      density_val =
        is_grain ? (_range_in[0] + rand_value * (_range_in[1] - _range_in[0]))
                  : (_range_out[0] + rand_value * (_range_out[1] - _range_out[0]));
    }
    
    elem_density[elem->id()] = density_val;

    std::vector<dof_id_type> dof_indices;
    dof_map.dof_indices(elem, dof_indices, var.number());

    //this explicit line assigns the computed values of density_val to the auxvariable supplied
    for (auto dof : dof_indices){
      sys.solution().set(dof, density_val);
    }
  }
  sys.solution().close();

  //read csv using density_i
  const std::vector<Real> data = _csv_total_fractions[0];
  const std::vector<Real> data_pore = _csv_total_fractions_pore[0];

  for (const auto & elem : _fe_problem.mesh().getMesh().active_element_ptr_range())
  {
    std::vector<dof_id_type> dof_indices_Y1;
    nl_dof_map.dof_indices(elem, dof_indices_Y1, var_Y1.number());

    if (dof_indices_Y1.empty())
      continue;

    // Elemental MicroID assigned earlier
    const Real density_value = elem_density[elem->id()];
    const int call_density = static_cast<int>(std::round(density_value));

    // Determine region type
    const bool is_pore = (call_density >= _range_pore[0] && call_density <= _range_pore[1]);
    const bool is_bulk = (call_density == _bulk_MicroID);
    const bool is_binder = (!is_pore && !is_bulk);

    // Assign fraction value
    Real predicted_fraction = 0.0;
    if (is_pore || is_bulk)
    {
      // For now, treat pores and bulk identically
      predicted_fraction = _bulk_RDX_fraction;
    }
    else
    {
      // Binder/interface: use CSV-based value
      const int idx = std::clamp(call_density, 0, static_cast<int>(data.size()) - 1);
      predicted_fraction = data[idx];
    }

    // Write the same fraction to all DOFs of this element
    for (auto dof : dof_indices_Y1)
      nl_sys.solution().set(dof, predicted_fraction);
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

//create a helper function to assign pore values

//second version
Real
PolycrystalDensityUO::assignPoreValue(unsigned int grain_id,
                                      const Elem * elem,
                                      bool is_grain,
                                      Real rand_value)
{
  // Static counter that persists during setup
  static std::unordered_map<unsigned int, unsigned int> pore_count_per_grain;

  bool is_pore = false;

  if (is_grain)
  {
    // deterministic pseudo-random number in [0,1) unique to elem+grain
    const Real rand_local =
        std::fmod(std::sin(elem->id() * 12.9898 + grain_id * 78.233) * 43758.5453, 1.0);

    // initialize counter for this grain the first time we see it
    if (pore_count_per_grain.find(grain_id) == pore_count_per_grain.end())
      pore_count_per_grain[grain_id] = 0;

    // only allow up to _n_pores per grain
    // pick the first _n_pores elements with smallest rand_local (< _n_pores / 1000 heuristic range)
    // or equivalently, use rand_local threshold but enforce upper bound via counter
    if (pore_count_per_grain[grain_id] < _n_pores)
    {
      // simple stochastic acceptance that tends to spread pores spatially
      const Real accept_prob = 0.002; // ~0.2% acceptance per candidate
      if (rand_local < accept_prob)
      {
        is_pore = true;
        pore_count_per_grain[grain_id]++;
      }
    }
  }

  // ---------- Assign value ----------
  Real returnvalue;
  if (is_pore)
  {
    // pick discrete pore type from _range_pore vector (e.g., 101, 102, 103)
    const unsigned int n_types = _range_pore.size();
    const unsigned int pore_index =
        static_cast<unsigned int>(std::floor(MooseRandom::rand() * n_types)) % n_types;
    returnvalue = static_cast<Real>(_range_pore[pore_index]);
  }
  else if (is_grain)
  {
    // normal grain material
    returnvalue = static_cast<Real>(_bulk_MicroID);
  }
  else
  {
    // binder / matrix
    returnvalue = _range_out[0] + rand_value * (_range_out[1] - _range_out[0]);
  }

  return returnvalue;
}


//Real
//PolycrystalDensityUO::assignPoreValue(unsigned int grain_id,
//                                      const Elem * elem,
//                                      bool is_grain,
//                                      Real rand_value)
//{
//
//  //cache pore id per grains
//  static std::unordered_map<unsigned int, std::vector<unsigned int>> pore_ids_per_grain;
//
//  //generate the rand pore set for each grain for the first time
//  if (pore_ids_per_grain.find(grain_id) == pore_ids_per_grain.end()){
//    std::vector<unsigned int> selected_ids;
//    selected_ids.reserve(_n_pores);
//
//    MooseRandom::seed(grain_id * 137 + 104729); //these numbers are randomly selected
//    std::unordered_set<unsigned int> unique_ids;
//    while (unique_ids.size() < _n_pores){
//      unsigned int candidate = static_cast<unsigned int>(MooseRandom::rand() * 1e6);
//      unique_ids.insert(candidate);
//    }
//    selected_ids.assign(unique_ids.begin(), unique_ids.end());
//    pore_ids_per_grain[grain_id] = selected_ids;
//  }
//
//  const unsigned int elem_id_hash = static_cast<unsigned int>(std::fmod(std::sin(elem->id() * 10 + grain_id * 70) * 30000, 1.) * 1e6);
//
//  bool is_pore = false;
//  const auto & selected_ids = pore_ids_per_grain[grain_id];
//  for (const auto & pid : selected_ids){
//    if (elem_id_hash == pid){
//      is_pore = true;
//      break;
//    }
//  }
//
//  Real returnvalue;
//  if (is_pore){
//    returnvalue = _range_pore[0] + rand_value * (_range_pore[1] - _range_pore[0]);
//  }
//  else if(is_grain){
//    returnvalue = _bulk_MicroID;
//  }
//  else{
//    returnvalue = _range_out[0] + rand_value * (_range_out[1] - _range_out[0]);
//  }
//  return returnvalue;
//}