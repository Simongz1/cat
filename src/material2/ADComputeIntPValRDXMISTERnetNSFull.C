#include "ADComputeIntPValRDXMISTERnetNSFull.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <fstream>
#include "Distribution.h"
#include "MooseRandom.h"

registerMooseObject("mlApp", ADComputeIntPValRDXMISTERnetNSFull);

InputParameters
ADComputeIntPValRDXMISTERnetNSFull::validParams()
{
  InputParameters params = Material::validParams();

  params.addClassDescription("compute MISTnetX predictions from interpolated values");

  //test: retrieve from name array
  params.addRequiredCoupledVar("v_components", "v_components");
  params.addRequiredCoupledVar("a_components", "a_components");

  params.addRequiredParam<Real>("thr_a", "acceleration threshold"); //deprecate this
  params.addRequiredParam<Real>("thr_v", "velocity threshold");

  //test: use gradient to compute activation

  params.addCoupledVar("density_i", "Coupled value");

  //CSV
  params.addRequiredParam<std::string>("csv_shock", "the name of the csv file with shock temperature values");
  params.addRequiredParam<std::string>("csv_react", "the name of the csv file with react temperature values");
  params.addRequiredParam<std::string>("csv_times", "the name of the csv file with times values");

  params.addRequiredParam<std::string>("csv_shock_pore", "csv_shock_pore");
  params.addRequiredParam<std::string>("csv_react_pore", "csv_react_pore");
  params.addRequiredParam<std::string>("csv_times_pore", "csv_times_pore");

  params.addRequiredParam<std::string>("csv_unreacted", "csv_unreacted");
  params.addRequiredParam<std::string>("csv_reacted", "csv_reacted");
  params.addRequiredCoupledVar("density_csv", "density_csv");

  //retrieve bulk ID
  params.addRequiredParam<unsigned int>("bulk_MicroID", "bulk_MicroID");
  params.addParam<unsigned int>("bulk_sensitivity", 20, "bulk_sensitivity");

  //retrieve range for pore IDS
  params.addRequiredParam<std::vector<unsigned int>>("range_pore", "range_pore");

  params.addRequiredParam<Real>("element_size", "element_size");
  params.addRequiredCoupledVar("tracking", "tracking");
  params.addParam<bool>("use_tabular_time", false, "use_tabular_time");
  params.addParam<bool>("use_distributions", false, "use_distributions");
  params.addParam<Real>("tau_react_scaling", 1., "tau_react_scaling");

  //add possible paramters for the us-up coefficients
  params.addParam<std::vector<Real>>("us_up_coeffs", {-0.286, 1.640, 1.249391, 5.575975}, "the coefficients for the polynomial fitting of the us-up relation from larger degree to smaller degree");
  return params;
}

ADComputeIntPValRDXMISTERnetNSFull::ADComputeIntPValRDXMISTERnetNSFull(
    const InputParameters & parameters)
  : Material(parameters),

  //rest
    _thr_a(getParam<Real>("thr_a")),
    _thr_v(getParam<Real>("thr_v")),
    _v_flag(declareADProperty<Real>("v_flag")),
    _v_flag_old(getMaterialPropertyOld<Real>("v_flag")),

    _temperature_mister_shock(declareProperty<Real>("temperature_mister_shock")),
    _temperature_mister_shock_old(getMaterialPropertyOld<Real>("temperature_mister_shock")),
    _temperature_mister_react(declareProperty<Real>("temperature_mister_react")),
    _temperature_mister_react_old(getMaterialPropertyOld<Real>("temperature_mister_react")),
    _density_i(coupledValue("density_i")),
    _called_up(declareADProperty<Real>("called_up")),
    _called_up_old(getMaterialPropertyOld<Real>("called_up")),
    _us(declareADProperty<Real>("us")),

    _csv_shock(getParam<std::string>("csv_shock")),
    _csv_react(getParam<std::string>("csv_react")),
    _csv_times(getParam<std::string>("csv_times")),

    //csv retrieval for pore
    _csv_shock_pore(getParam<std::string>("csv_shock_pore")),
    _csv_react_pore(getParam<std::string>("csv_react_pore")),
    _csv_times_pore(getParam<std::string>("csv_times_pore")),

    //declare time
    _time_react(declareADProperty<Real>("time_react")),
    _time_react_old(getMaterialPropertyOld<Real>("time_react")),

    _csv_unreacted(getParam<std::string>("csv_unreacted")),
    _csv_reacted(getParam<std::string>("csv_reacted")),
    _density(declareADProperty<Real>("density")),
    _density_csv(coupledValue("density_csv")),

    //bulk grains assignment
    _bulk_MicroID(getParam<unsigned int>("bulk_MicroID")),
    _bulk_sensitivity(getParam<unsigned int>("bulk_sensitivity")),
    _range_pore(getParam<std::vector<unsigned int>>("range_pore")),

    _rate_tracking(declareADProperty<Real>("rate_tracking")),
    _h(getParam<Real>("element_size")),
    _tracking(adCoupledValue("tracking")),
    _use_tabular_time(getParam<bool>("use_tabular_time")),
    _use_distributions(getParam<bool>("use_distributions")),

    //placeholder for distribution call
    _distribution_lower(nullptr),
    _distribution_upper(nullptr),
    _tau_react_scaling(getParam<Real>("tau_react_scaling")),
    //get parameters for us-up
    _coeffs(getParam<std::vector<Real>>("us_up_coeffs"))
{ 
  /////////////////////////////////////////////////////////

  const unsigned int n_v = coupledComponents("v_components");
  _v.reserve(n_v);
  for (unsigned int i = 0; i < n_v; ++i){
    _v.push_back(&adCoupledValue("v_components", i));
  }
    
  const unsigned int n_a = coupledComponents("a_components");
  _a.reserve(n_a);
  for (unsigned int i = 0; i < n_a; ++i){
    _a.push_back(&adCoupledValue("a_components", i));
  }

  //make the total CSV local
  const auto csv_total_shock = readCSV(_csv_shock);
  const auto csv_total_react = readCSV(_csv_react);
  
  if (_use_tabular_time){
    _csv_total_times = readCSV(_csv_times);
  }

  //read all the csv structures that are needed
  const auto csv_total_shock_pore = readCSV(_csv_shock_pore);
  const auto csv_total_react_pore = readCSV(_csv_react_pore);
  const auto csv_total_times_pore = readCSV(_csv_times_pore);
  const auto csv_total_pu = readCSV(_csv_unreacted);
  const auto csv_total_pr = readCSV(_csv_reacted);

  //use the csv function to form the velocity vectors
  _up_values = formVelocityVector(csv_total_shock);
  _up_values_pore = formVelocityVector(csv_total_shock_pore);

  //retrieve nanoPBX csv data using the new helper function
  _temperature_values_shock = formData(csv_total_shock);
  _temperature_values_react = formData(csv_total_react);

  if (_use_tabular_time){
    _time_values = formData(_csv_total_times);

  }

  //retrieve csv data for pore
  _temperature_values_shock_pore = formData(csv_total_shock_pore);
  _temperature_values_react_pore = formData(csv_total_react_pore);
  _time_values_pore = formData(csv_total_times_pore);
}

void 
ADComputeIntPValRDXMISTERnetNSFull::initQpStatefulProperties()
{
  _v_flag[_qp] = ADReal(0.0);
  _stored_shock = _temperature_mister_shock_old[_qp];
  _stored_react = _temperature_mister_react_old[_qp];
  _time_react[_qp]  = ADReal(_time_react_old[_qp]);
  _called_up[_qp] = _called_up_old[_qp];
  _rate_tracking[_qp] = 0.;
}

void
ADComputeIntPValRDXMISTERnetNSFull::computeQpProperties()
{
  //store density in material property
  _density[_qp] = ADReal(_density_csv[_qp]);
  //KEEP THIS ORDER

  //ORDER 1
  //keep stateful value 1
  _v_flag[_qp] = ADReal(_v_flag_old[_qp]);
  _temperature_mister_shock[_qp] = _temperature_mister_shock_old[_qp];
  _temperature_mister_react[_qp] = _temperature_mister_react_old[_qp];
  _time_react[_qp] = ADReal(_time_react_old[_qp]);
  _called_up[_qp] = _called_up_old[_qp];

  //ORDER 2
  //activate shock heat: call when velocity is bigger than a value 1

  //build velocity and acceleration vectors
  std::vector<ADReal> v_vect(_v.size());
  std::vector<ADReal> a_vect(_a.size());

  for (unsigned int i = 0; i < _v.size(); ++i){
    v_vect[i] = (*_v[i])[_qp];
  }
  for (unsigned int j = 0; j < _a.size(); ++j){
    a_vect[j] = (*_a[j])[_qp];
  }

  //define inline for L2norm
  auto L2norm = [](const std::vector<ADReal> & vect) -> ADReal
  {
    ADReal sum = 0;
    for (int i = 0; i < vect.size(); ++i){
      sum += (vect[i] * vect[i]);
    }
    return MetaPhysicL::sqrt(sum);
  };

  auto time_equil = [this](const ADReal & up) -> ADReal
  {
    //define condition for deflagrating or quench
    if (_temperature_mister_react[_qp] > 1000){
      ADReal time_def = 8.046150413374 * MetaPhysicL::exp(-0.963962928923 * up);
      return time_def;
    }
    else{
      ADReal time_quench = 0.000082188048 * MetaPhysicL::exp(3.430175735695 * up);
      return time_quench;
    }
  };

  ADReal v_norm;
  ADReal a_norm;
  v_norm = L2norm(v_vect);
  a_norm = L2norm(a_vect);

  //retrieve polynomial function
  //in the future, this should read the parameters from the input
  //_us[_qp] = -0.286 * MetaPhysicL::pow(v_norm, 3.) + 
  //            1.640 * MetaPhysicL::pow(v_norm, 2.) -
  //            1.249391 * MetaPhysicL::pow(v_norm, 1.) +
  //            5.575975;
  
  //use the helper function
  _us[_qp] = computeUs(v_norm, _coeffs);
  
  //declare rate for shock tracking
  _rate_tracking[_qp] = (v_norm > _thr_v) ? _us[_qp] / (3 * _h) : ADReal(0.); //this is the time it takes for the shock to travel 3 elements
  //form call condition by evaluating when v is equal to a

  //LOGIC:
  //1. Expect a to be greater than a value, around 1
  //2. When norm(v) = norm(a) we call
  //3. Call is still stateful

  bool _call_condition = false;

  //to get the old behavior, use up/us and thr_a to be smaller than a value

  if (_qp==0){
    if(_v_flag[_qp]==0.){
      if (v_norm >= _thr_v && _tracking[_qp] >= ADReal(0.99)){
        _call_condition = true;
      }
    }
  }

  if (_call_condition){ //we need V and A constraints to make sure the call happens at the actual shock velocity
    _v_flag[_qp] = 1.0; //set flag to one, call in misternet material as condition for ComputeQpProperties
    //get temperatures
    const int id_call = static_cast<int>(std::round(_density_i[_qp]));

    //here we need to branch if we are using bulk or MicroID values

    //declare values for binder
    Real pred_shock, pred_react, pred_time;

    //define cases
    //cases for when loaded_micro = false
    bool is_pore = false;
    bool is_bulk = false;
    bool is_binder = false;

    //standard distribution case
    is_pore = (id_call >= _range_pore[0] && id_call <= _range_pore[1]);
    is_bulk = (id_call == _bulk_MicroID);
    is_binder = (!is_pore && !is_bulk);
    
    //define call velocity explicitly
    const Real call_up = std::clamp(L2norm(v_vect).value(), (0.0), (4.89));

    //retrieve data based on grain or binder
    if (is_binder){ //this is the usual loop
      //assign hand coded values as placeholders
      //call the pair of values only once
      const std::pair<Real, Real> temps = getTemperatures(call_up, id_call, "binder");
      pred_shock = temps.first;
      pred_react = temps.second;
      if (_use_tabular_time){
        pred_time  = getTimes(call_up, id_call, "binder");
      }
    }
    else if (is_pore){ //this is the loop for grains (pore + bulk)
      //manually clamp the id
      int id_call_pore;
      if (id_call >= _range_pore.back()){
        id_call_pore = _range_pore.back();
      }
      if (id_call <= _range_pore.front()){
        id_call_pore = _range_pore.front();
      }
      //call temperature only once and store both values
      const std::pair<Real, Real> temps_pore = getTemperatures(call_up, id_call_pore, "pore");
      
      pred_shock = temps_pore.first;
      pred_react = temps_pore.second;
      pred_time  = getTimes(call_up, id_call_pore, "pore");
    }
    else{
      //make consistent here too
      const std::pair<Real, Real> temps_bulk = getTemperatures(call_up, id_call, "bulk");
      pred_shock = temps_bulk.first;
      pred_react = temps_bulk.second;

      if (_use_tabular_time){
        pred_time  = getTimes(call_up, id_call, "bulk");
      }
      //
    }
  
    //store in material property
    _temperature_mister_shock[_qp] = pred_shock;
    _temperature_mister_react[_qp] = pred_react;

    if (_use_tabular_time){
      _time_react[_qp] = pred_time;
    }else if (_use_distributions){
      _time_react[_qp] = getDistributionTime(call_up, _temperature_mister_react[_qp]);
      _time_react[_qp] = max(_dt, _time_react[_qp]);
    }else{
      _time_react[_qp] = time_equil(ADReal(call_up));
      _time_react[_qp] = max(_dt, min(2.0, _time_react[_qp]));
    }

    ///scale tau
    _time_react[_qp] *= _tau_react_scaling;

    //get the temeprature at the initial qp
    _stored_shock = pred_shock;
    _stored_react = pred_react;

    if (_use_tabular_time){
      _stored_time  = pred_time;
    }
    
    _called_up[_qp] = L2norm(v_vect);
  }

  //ORDER 3
  if (_qp > 0.){
    _v_flag[_qp] = _v_flag[0];
    _temperature_mister_shock[_qp] = _temperature_mister_shock[0];
    _temperature_mister_react[_qp] = _temperature_mister_react[0];
    _time_react[_qp] = _time_react[0];
    _called_up[_qp] = _called_up[0];
  }  
}

//interpolate between values
//use pointers to avoid memory issues
//make interpolation term wise 
Real
ADComputeIntPValRDXMISTERnetNSFull::interpolation(const std::vector<Real> & A, const std::vector<Real> & B, const Real t, const unsigned int index){
  //directly use the index to perform interpolation
  //std::vector<Real> res;
  //res.reserve(A.size());
  //for (size_t i = 0; i < A.size(); ++i){
  //  res.push_back((1. - t) * A[i] + t * B[i]);
  //}
  Real value;
  value = (1. - t) * A[index] + t * B[index];
  return value;
}

//helper function to read CSV file
std::vector<std::vector<Real>>
ADComputeIntPValRDXMISTERnetNSFull::readCSV(const std::string csv_file_name){
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

//helper function to retrieve temperatures
std::pair<Real, Real>
ADComputeIntPValRDXMISTERnetNSFull::getTemperatures(const Real up, const int id, const std::string phase){

  //forward declaration
  Real lower_bound, upper_bound, interval_number;

  //change these to null pointers
  const std::vector<Real> * lower_temps_shock = nullptr;
  const std::vector<Real> * upper_temps_shock = nullptr;
  const std::vector<Real> * lower_temps_react = nullptr;
  const std::vector<Real> * upper_temps_react = nullptr;

  //std::vector<Real> lower_temps_shock, upper_temps_shock;
  //std::vector<Real> lower_temps_react, upper_temps_react;

  //to use the same interval, call it here
  _interval = computeInterval(up, _up_values);
  _interval_pore = computeInterval(up, _up_values_pore);

  //here we need to branch based on the phase
  if (phase == "binder"){
    //direct assignment of the shock temperatures
    //dereference the pointers here to assign value
    lower_temps_shock = &_temperature_values_shock[_interval];
    upper_temps_shock = &_temperature_values_shock[_interval + 1];

    //direct assignment of the reaction temperatures
    lower_temps_react = &_temperature_values_react[_interval];
    upper_temps_react = &_temperature_values_react[_interval + 1];

    //and also assign bounds for interpolation
    lower_bound = _up_values[_interval];
    upper_bound = _up_values[_interval + 1];
  }
  else if (phase == "pore"){
    //direct assignment of shock temperatures
    lower_temps_shock = &_temperature_values_shock_pore[_interval_pore];
    upper_temps_shock = &_temperature_values_shock_pore[_interval_pore + 1];

    //direct assignment of reaction temperatures
    lower_temps_react = &_temperature_values_react_pore[_interval_pore];
    upper_temps_react = &_temperature_values_react_pore[_interval_pore + 1];

    //assign bounds
    lower_bound = _up_values_pore[_interval_pore];
    upper_bound = _up_values_pore[_interval_pore + 1];
  }
  else if (phase == "bulk"){
    //direct assignment of shock temperatures
    lower_temps_shock = &_temperature_values_shock[_interval];
    upper_temps_shock = &_temperature_values_shock[_interval + 1];

    //direct assingmnet of reaction temperatures
    lower_temps_react = &_temperature_values_react[_interval];
    upper_temps_react = &_temperature_values_react[_interval + 1];

    //assign bounds
    lower_bound = _up_values[_interval];
    upper_bound = _up_values[_interval + 1];
  }

  //this assigns the ratio
  Real t = (up - lower_bound) / (upper_bound - lower_bound);

  //test: use this same t to interpolate time
  _ratio = t;

  Real temp_shock;
  Real temp_react;

  if (phase == "binder"){ //traditional branch
    temp_shock = interpolation(*lower_temps_shock, *upper_temps_shock, t, id);
    temp_react = interpolation(*lower_temps_react, *upper_temps_react, t, id);
  }

  if (phase == "pore"){
    //the pore microstructures are 101, 102, 103, so we need to subtract to access the data
    int id_pore = static_cast<int>(std::round(id - _range_pore.front()));

    id_pore = std::max(0, std::min(id_pore, static_cast<int>(lower_temps_shock->size()) - 1));
    
    //now access pore data
    temp_shock = interpolation(*lower_temps_shock, *upper_temps_shock, t, id_pore);
    temp_react = interpolation(*lower_temps_react, *upper_temps_react, t, id_pore);
  }
  if (phase == "bulk"){
    //bulk will preserve the nanoPBX data, but will be called all at an insensitive microstructure
    int id_bulk = _bulk_sensitivity; //the most insensitive, this can be tuned until bulk data is available
    
    //now access bulk data
    temp_shock = interpolation(*lower_temps_shock, *upper_temps_shock, t, id_bulk);
    temp_react = interpolation(*lower_temps_react, *upper_temps_react, t, id_bulk);
  }

  //make the return type consistent
  return std::make_pair(temp_shock, temp_react);
}

Real
ADComputeIntPValRDXMISTERnetNSFull::getTimes(const Real up, const int id, const std::string phase){

  //create times lower and upper
  const std::vector<Real> * lower_times;
  const std::vector<Real> * upper_times;

  //define lower and upper times
  //here we account for the phase of the element

  if (phase == "binder"){
    //traditional branch
    lower_times = &_time_values[_interval];
    upper_times = &_time_values[_interval + 1];
  }
  if (phase == "pore"){
    //here we use the pore time values
    lower_times = &_time_values_pore[_interval];
    upper_times = &_time_values_pore[_interval + 1];
  }
  if (phase == "bulk"){
    //preserve base nanoPBX
    lower_times = &_time_values[_interval];
    upper_times = &_time_values[_interval + 1];
  }

  //interpolate based on previously computed interval and balance number

  //return based on microstructures
  Real time_value;

  if (phase == "binder"){
    time_value = interpolation(*lower_times, *upper_times, _ratio, id);
  }

  if (phase == "pore"){
    const int id_pore = id - 101;
    time_value = interpolation(*lower_times, *upper_times, _ratio, id_pore);
  }

  if (phase == "bulk"){
    int id_bulk = _bulk_sensitivity;
    time_value = interpolation(*lower_times, *upper_times, _ratio, id_bulk);
  }

  //return whatever value went through
  return time_value;
}

//function for retrieving distribution data
Real
ADComputeIntPValRDXMISTERnetNSFull::getDistributionTime(const Real up, const Real predicted_temp){
  //with the known call up, call the closest lower and upper distributions
  const Real lower_up = std::floor(10. * up);
  const Real upper_up = std::ceil(10. * up);

  //scale up
  const Real scaled_up = 10. * up;

  //determine whether it is a deflagrated or quenched microstructure
  const bool deflagrated = predicted_temp > 1000. ? true : false;

  //form the names
  const std::string appendix = deflagrated ? "_reacted" : "_unreacted";
  const DistributionName lower_name = std::to_string(static_cast<int>(lower_up)) + appendix;
  const DistributionName upper_name = std::to_string(static_cast<int>(upper_up)) + appendix;

  //retrieve the lower and upper distributions
  const Distribution & lower_distribution = getDistributionByName(lower_name);
  const Distribution & upper_distribution = getDistributionByName(upper_name);

  //seed number to call both distributions
  const Real seed = MooseRandom::rand();

  //evaluate the distribution with a randomly generated number
  const Real lower_time = lower_distribution.quantile(seed);
  const Real upper_time = upper_distribution.quantile(seed);

  //interpolate based on distance
  const Real dist = (scaled_up - lower_up) / (upper_up - lower_up);
  const Real interpolated_value = lower_time + (dist) * (upper_time - lower_time);

  return interpolated_value;
}

//we can define an explicit function that computes the interval of u_p where we are

Real
ADComputeIntPValRDXMISTERnetNSFull::computeInterval(const Real up, const std::vector<Real> & up_values){
  //define value to retunr
  //dereference the pointer for the vector of up values
  const std::vector<Real> & up_vals = up_values;

  //iterate
  Real interval_number = 0.;
  for (unsigned int i = 1; i < up_vals.size(); ++i){
    //find the i values where up_call is between the lower and the upper bounds of the interval
    if (up_vals[i] > up){
      interval_number = i - 1;
      break;
    }
  }
  return interval_number;
}

//define a function that forms the vectors neeeded for interpolation
//this takes vector of vectors and returns the member vector

std::vector<std::vector<Real>>
ADComputeIntPValRDXMISTERnetNSFull::formData(const std::vector<std::vector<Real>> & csv){
  //form output
  std::vector<std::vector<Real>> data;

  //iterate
  for (auto &row : csv){
    if (row.size() < 2){
      mooseError("need bigger CSV");
    }
  
    data.emplace_back(row.begin() + 1, row.end());
  }
  return data;
}

//also generate a function to form velocity vectors for intepolation
std::vector<Real>
ADComputeIntPValRDXMISTERnetNSFull::formVelocityVector(const std::vector<std::vector<Real>> & csv){
  //form output
  std::vector<Real> velocity;
  for (auto & row : csv){
    //no need to check for size here
    velocity.emplace_back(row[0]);
  }
  return velocity;
}

//define helper function to evaluate us-up relation
ADReal
ADComputeIntPValRDXMISTERnetNSFull::computeUs(const ADReal & up, const std::vector<Real> & coeffs){
  //form the polynomial from the coefficients
  ADReal us = ADReal(0.);
  for (unsigned int i = 0; i < coeffs.size(); ++i){
    us += coeffs[i] * MetaPhysicL::pow(up, static_cast<int>(coeffs.size() - 1 - i));
  }
  return us;
}