#include "PolycrystalDensityUO.h"
#include "Function.h"

registerMooseObject("mlApp", PolycrystalDensityUO);

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
  params.addParam<Real>("pore_probability", 0.01, "pore_probability"); //defaul value
  params.addParam<bool>("euler_angles", true, "euler_angles");

  //for loaded custom microstructure
  params.addParam<bool>("use_loaded_microstructure", false, "use_loaded_microstructure");
  //name to be used to directly retrieve the function
  params.addParam<FunctionName>("loaded_microstructure", "loaded_microstructure", "loaded_microstructure");
  params.addParam<Real>("pore_limit", 0.25, "pore_limit");
  params.addParam<Real>("bulk_limit", 1.77, "bulk_limit");
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
    _pore_RDX_fraction(getParam<Real>("pore_RDX_fraction")),
    _pore_probability(getParam<Real>("pore_probability")),
    _euler_angles(getParam<bool>("euler_angles")),

    //for loaded microstructure
    _use_loaded_microstructure(getParam<bool>("use_loaded_microstructure")),
    _loaded_microstructure_name(getParam<FunctionName>("loaded_microstructure")),
    _pore_limit(getParam<Real>("pore_limit")),
    _bulk_limit(getParam<Real>("bulk_limit"))
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

  //allocate centers
  _centers.clear();
  _centers.reserve(_num_grains);

  //allocate radii
  _radii.clear();
  _radii.reserve(_num_grains);

  //allocate grainID
  _grainID.clear();
  _grainID.reserve(_num_grains);

  //get the function that contains the data
  const Function & loaded_function = getFunction(_loaded_microstructure_name);

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
  ////////////////////////
  if (_tid == 0){
     mooseInfo("Generated ", _num_grains, " Voronoi centers in PolycrystalDensityUO. Next step is to assign defects inside grains and nanoPBXs at the interfaces");
  }
  ///////////////////////

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

          if (r < _pore_probability){ //this threshold can be changed to generate more pore sites
            pores_per_grain[gid].insert(elem);
            count++;
          }
        }
      }
    }
  }

  //here we need to generate cracks in a similar way as we generate the pores
  
  //here starts the variable assignment
  auto & sys = _fe_problem.getAuxiliarySystem();
  auto & nl_sys = _fe_problem.getNonlinearSystem(0);
  auto & var = sys.getVariable(_tid, "density_i");
  auto & var_grainID = sys.getVariable(_tid, "grainID");
  auto & var_Y1 = nl_sys.getVariable(_tid, "Y1");
  const DofMap & dof_map = sys.system().get_dof_map();
  const DofMap & nl_dof_map = nl_sys.system().get_dof_map();

  //loaded microstructure variable
  //auto & loaded_microstructure = sys.getVariable(_tid, _loaded_microstructure_name);

  auto & var_euler1 = sys.getVariable(_tid, "euler1");
  auto & var_euler2 = sys.getVariable(_tid, "euler2");
  auto & var_euler3 = sys.getVariable(_tid, "euler3");

  std::unordered_map<dof_id_type, Real> elem_euler1;
  std::unordered_map<dof_id_type, Real> elem_euler2;
  std::unordered_map<dof_id_type, Real> elem_euler3;

  //generate euler angles vectors
  std::vector<Real> grain_euler1(_num_grains), grain_euler2(_num_grains), grain_euler3(_num_grains);

  if (_euler_angles){
    MooseRandom::seed(1234);

    for (unsigned int g = 0; g < _num_grains; ++g){
      grain_euler1[g] = MooseRandom::rand();
      grain_euler2[g] = MooseRandom::rand();
      grain_euler3[g] = MooseRandom::rand();
    }

    if (_tid == 0){
      mooseInfo("Generated random euler fractions");
    }
  }

  //unordered map for density
  std::unordered_map<dof_id_type, Real> elem_density;

  //unordered map for grainID
  std::unordered_map<dof_id_type, Real> elem_grainID;

  //write initial values to the variable field
  for (const auto & elem : _fe_problem.mesh().getMesh().active_element_ptr_range())
  {
    Point centroid = elem->vertex_average();

    Real min_dist = std::numeric_limits<Real>::max();
    Real second_min_dist = std::numeric_limits<Real>::max();
    unsigned int nearest = 0;
    
    //this finds the nearest center at each element
    //for (unsigned int i = 0; i < _centers.size(); ++i)
    //{
    //  const Real d = (centroid - _centers[i]).norm();
    //  if (d < min_dist)
    //  {
    //    second_min_dist = min_dist;
    //    min_dist = d;
    //    nearest = i;
//
    //    //use the nearest center to assign the grainID
    //    _grainID.push_back(i + 1); //this has an arbitrary reference at 1, 0 will be left for binder
//
    //    //now we need to make sure to assign grain ID only from 1 to _n_grains, 0 will be binder
    //  }
    //  else if (d < second_min_dist){
    //    second_min_dist = d;
    //  }
    //}

    //generalized definition of three cases
    bool is_pore = false;
    bool is_bulk = false;
    bool is_binder = false;
    bool is_grain = false;
    bool is_boundary = false;
    bool is_far = false;

    const Real boundary_gap = (second_min_dist - min_dist);
    
    //
    const bool in_target =
        std::find(_target_grains.begin(), _target_grains.end(), nearest + 1) != _target_grains.end();

    const Real rand_value = MooseRandom::rand();

    //create branch for the case when _bulk_grains = true
    //this assigns an artificial/placeholder MicroID to all elements inside grains
    Real density_val;

    //create local variable for grainID to use to assign values later
    Real grainID;

    //declare real assignments
    Real euler1 = 0.;
    Real euler2 = 0.;
    Real euler3 = 0.;

    Real loaded_density_value = 0.0;
    if (_bulk_grains){ //use_loaded_microstructure needs this
      //do the branch here
      if (_use_loaded_microstructure){
        //read the function and evaluate at the centroid of each element
        //const Function & loaded_function = getFunction(_loaded_microstructure_name);
        loaded_density_value = loaded_function.value(0.0, centroid);

        //define static cases
        is_pore = loaded_density_value <= _pore_limit;
        is_bulk = loaded_density_value >= _bulk_limit;
        is_binder = !is_pore && !is_bulk;
      }else{
        //standard approach
        is_far = min_dist > _radii[nearest];
        is_boundary = boundary_gap < _matrix_thickness;
        is_grain = !is_boundary && !is_far;
        
        //old assignment
        is_pore = pores_per_grain.count(nearest) && pores_per_grain[nearest].count(elem);
        is_bulk = (is_grain && !is_pore);
        is_binder = (!is_pore && !is_bulk);
      }

      //now use the defined cases
      if (is_pore){
        //standard pore assignment
        const unsigned int n_types = _range_pore.size();
        const unsigned int idx = static_cast<unsigned int>(std::floor(MooseRandom::rand() * n_types)) % n_types;
        
        //density
        density_val = static_cast<Real>(_range_pore[idx]);
        grainID = static_cast<Real>(nearest + 1);
      }else if (is_bulk){
        density_val = static_cast<Real>(_bulk_MicroID);
        grainID = static_cast<Real>(nearest + 1);

        //assign euler
        euler1 = grain_euler1[nearest];
        euler2 = grain_euler2[nearest];
        euler3 = grain_euler3[nearest];
      }
      
      //else if (is_grain){
      //  density_val = static_cast<Real>(_bulk_MicroID);
      //  grainID = static_cast<Real>(nearest + 1);

      //  euler1 = grain_euler1[nearest];
      //  euler2 = grain_euler2[nearest];
      //  euler3 = grain_euler3[nearest];
      //}

      else{
        density_val = _range_out[0] + rand_value * (_range_out[1] - _range_out[0]);
        grainID = 0; //this corresponds to binder
      }
    }else{
      //older approach
      is_far = min_dist > _radii[nearest];
      is_boundary = boundary_gap < _matrix_thickness;
      is_grain = !is_boundary && !is_far;

      //assign value
      density_val =
        is_grain ? (_range_in[0] + rand_value * (_range_in[1] - _range_in[0]))
                  : (_range_out[0] + rand_value * (_range_out[1] - _range_out[0]));
    }
    
    //here the density value is assigned to the unordered map
    elem_density[elem->id()] = density_val;
    elem_grainID[elem->id()] = grainID;

    elem_euler1[elem->id()] = euler1;
    elem_euler2[elem->id()] = euler2;
    elem_euler3[elem->id()] = euler3;

    //this is specific for density, replicate for grainID
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_grainID;

    //generate dof map if euler angles is set
    if (_euler_angles){
      std::vector<dof_id_type> dof_indices_euler1;
      std::vector<dof_id_type> dof_indices_euler2;
      std::vector<dof_id_type> dof_indices_euler3;

      dof_map.dof_indices(elem, dof_indices_euler1, var_euler1.number());
      dof_map.dof_indices(elem, dof_indices_euler2, var_euler2.number());
      dof_map.dof_indices(elem, dof_indices_euler3, var_euler3.number());

      //write into variables

      for (auto dof : dof_indices_euler1){
        sys.solution().set(dof, euler1);
      }
      for (auto dof : dof_indices_euler2){
        sys.solution().set(dof, euler2);
      }
      for (auto dof : dof_indices_euler3){
        sys.solution().set(dof, euler3);
      }
    }

    dof_map.dof_indices(elem, dof_indices, var.number());
    dof_map.dof_indices(elem, dof_indices_grainID, var_grainID.number());

    //this explicit line assigns the computed values of density_val to the auxvariable supplied
    for (auto dof : dof_indices){
      sys.solution().set(dof, density_val);
    }
    for (auto dof : dof_indices_grainID){
      sys.solution().set(dof, grainID);
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

    //elemental MicroID assigned earlier
    const Real density_value = elem_density[elem->id()];
    const int call_density = static_cast<int>(std::round(density_value));

    //region type conditions
    const bool is_pore = (call_density >= _range_pore[0] && call_density <= _range_pore[1]);
    const bool is_bulk = (call_density == _bulk_MicroID);
    const bool is_binder = (!is_pore && !is_bulk);

    //fraction value based on region conditions
    Real predicted_fraction = 0.0;
    if (is_pore || is_bulk)
    {
      //here we assume that pores are 100%RDX MASS FRACTION
      predicted_fraction = _bulk_RDX_fraction;
    }
    else
    {
      //binder and nanoPBXs
      const int idx = std::clamp(call_density, 0, static_cast<int>(data.size()) - 1);
      predicted_fraction = data[idx];
    }

    //write fraction
    for (auto dof : dof_indices_Y1)
      nl_sys.solution().set(dof, predicted_fraction);
  }

  nl_sys.solution().close();
}
///kept empty on purpose
void
PolycrystalDensityUO::execute()
{}
void
PolycrystalDensityUO::initialize()
{}
void 
PolycrystalDensityUO::finalize()
{}

//standard function to read CSV data
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
//this should also account for the assignment when loaded is used
Real
PolycrystalDensityUO::assignPoreValue(unsigned int grain_id,
                                      const Elem * elem,
                                      bool is_grain,
                                      Real rand_value)
{
  //define pore counter per grain
  static std::unordered_map<unsigned int, unsigned int> pore_count_per_grain;

  //initially assume that the current element is not a pore
  bool is_pore = false;

  if (is_grain)
  {
    //generate a random based on the element ID, this should be tuneable from the input file
    const Real rand_local =
        std::fmod(std::sin(elem->id() * 12.9898 + grain_id * 78.233) * 43758.5453, 1.0);

    //this means that the element is the last on the appended list => first time seeing it
    if (pore_count_per_grain.find(grain_id) == pore_count_per_grain.end())
      pore_count_per_grain[grain_id] = 0;

    //only allow up to _n_pores per grain
    //here we use the _n_pores parameter
    if (pore_count_per_grain[grain_id] < _n_pores)
    {
      //here we also apply a probabilistic pore acceptance fraction
      const Real accept_prob = 0.002;  //this one also should come from the input file
      if (rand_local < accept_prob)
      {
        is_pore = true;
        pore_count_per_grain[grain_id]++;
      }
    }
  }

  //value assignments after generation
  Real returnvalue;
  if (is_pore)
  {
    const unsigned int n_types = _range_pore.size();
    const unsigned int pore_index =
        static_cast<unsigned int>(std::floor(MooseRandom::rand() * n_types)) % n_types;
    returnvalue = static_cast<Real>(_range_pore[pore_index]);
  }
  else if (is_grain)
  {
    returnvalue = static_cast<Real>(_bulk_MicroID);
  }
  else
  {
    returnvalue = _range_out[0] + rand_value * (_range_out[1] - _range_out[0]);
  }

  return returnvalue;
}
////////////////////////////////