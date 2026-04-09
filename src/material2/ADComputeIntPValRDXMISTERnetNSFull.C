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

  params.addClassDescription("Standard compute Mie Gruneisen Pressure with JWL pressure for reacted material. Also computes artificial viscosity contribution");

  //test: retrieve from name array
  params.addRequiredCoupledVar("v_components", "v_components");
  params.addRequiredCoupledVar("a_components", "a_components");

  params.addRequiredParam<Real>("thr_a", "acceleration threshold");
  params.addRequiredParam<Real>("thr_v", "velocity threshold");

  //test: use gradient to compute activation

  params.addCoupledVar("density_i", "Coupled value");

  //CSV
  params.addRequiredParam<std::string>("csv_shock", "csv_shock");
  params.addRequiredParam<std::string>("csv_react", "csv_react");
  params.addRequiredParam<std::string>("csv_times", "csv_times");

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
    _tau_react_scaling(getParam<Real>("tau_react_scaling"))
{ 
  /////////////////////////////////////////////////////////
  const unsigned int n_v = coupledComponents("v_components");
  _v.reserve(n_v);
  for (unsigned int i = 0; i < n_v; ++i)
    _v.push_back(&adCoupledValue("v_components", i));

  const unsigned int n_a = coupledComponents("a_components");
  _a.reserve(n_a);
  for (unsigned int i = 0; i < n_a; ++i)
    _a.push_back(&adCoupledValue("a_components", i));
  
  /////////////////////////////////////////////////////////

  /////%%%%%%%%%%%%%%%%%%%%%%%%%/////////////////////////%%%%%%%//////
  //here I cache the table only once
  _csv_total_shock = readCSV(_csv_shock);
  _csv_total_react = readCSV(_csv_react);

  if (_use_tabular_time){
    _csv_total_times = readCSV(_csv_times);
  }

  _csv_total_shock_pore = readCSV(_csv_shock_pore);
  _csv_total_react_pore = readCSV(_csv_react_pore);
  _csv_total_times_pore = readCSV(_csv_times_pore);
  //
  _csv_total_pu = readCSV(_csv_unreacted);
  _csv_total_pr = readCSV(_csv_reacted);

  //retrieve nanoPBX csv data

  for (auto &row : _csv_total_shock){
    if (row.size() < 2){
      mooseError("need bigger CSV");
    }
    _up_values.push_back(row[0]);
    std::vector<Real> temps(row.begin() + 1, row.end());
    _temperature_values_shock.push_back(std::move(temps));
  }

  for (auto &row : _csv_total_react){
    if (row.size() < 2){
      mooseError("need bigger CSV");
    }
    std::vector<Real> temps(row.begin() + 1, row.end());
    _temperature_values_react.push_back(std::move(temps));
  }

  if (_use_tabular_time){
    for (auto &row : _csv_total_times){
      if (row.size() < 2){
        mooseError("need bigger CSV");
      }
      std::vector<Real> times(row.begin() + 1, row.end());
      _time_values.push_back(std::move(times));
    }
  }
  
  //retrieve csv data for pore
  for (auto &row : _csv_total_shock_pore){
    if (row.size() < 2){
      mooseError("need bigger CSV");
    }
    _up_values_pore.push_back(row[0]);
    std::vector<Real> temps(row.begin() + 1, row.end());
    _temperature_values_shock_pore.push_back(std::move(temps));
  }

  for (auto &row : _csv_total_react_pore){
    if (row.size() < 2){
      mooseError("need bigger CSV");
    }
    std::vector<Real> temps(row.begin() + 1, row.end());
    _temperature_values_react_pore.push_back(std::move(temps));
  }

  for (auto &row : _csv_total_times_pore){
    if (row.size() < 2){
      mooseError("need bigger CSV");
    }
    std::vector<Real> times(row.begin() + 1, row.end());
    _time_values_pore.push_back(std::move(times));
  }

  ////////////////////////////////////////////////////////////////
  //the csv structure is J, P_u, P_r
  //this gets the pressures as a list on a vector

  ////////////////////////////////////////////////////////////////
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
  auto L2norm = [](const std::vector<ADReal> &vect) -> ADReal
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

  ADReal condition_v;
  ADReal condition_a;
  condition_v = L2norm(v_vect);
  condition_a = L2norm(a_vect);

  _us[_qp] = -0.286 * MetaPhysicL::pow(condition_v, 3.) + 
              1.640 * MetaPhysicL::pow(condition_v, 2.) -
              1.249391 * MetaPhysicL::pow(condition_v, 1.) +
              5.575975;
  
  //declare rate for shock tracking
  _rate_tracking[_qp] = (condition_v > _thr_v) ? _us[_qp] / (3 * _h) : ADReal(0.); //this is the time it takes for the shock to travel 3 elements
  //form call condition by evaluating when v is equal to a

  //LOGIC:
  //1. Expect a to be greater than a value, around 1
  //2. When norm(v) = norm(a) we call
  //3. Call is still stateful

  bool _call_condition = false;

  //to get the old behavior, use up/us and thr_a to be smaller than a value

  if (_qp==0){
    if(_v_flag[_qp]==0.){
      if (condition_v >= _thr_v && _tracking[_qp] >= ADReal(0.99)){
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
      pred_shock = getTemperatures(call_up, id_call, "binder")[0];
      pred_react = getTemperatures(call_up, id_call, "binder")[1];
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

      pred_shock = getTemperatures(call_up, id_call_pore, "pore")[0];
      pred_react = getTemperatures(call_up, id_call_pore, "pore")[1];
      pred_time  = getTimes(call_up, id_call_pore, "pore");
    }
    else{
      pred_shock = getTemperatures(call_up, id_call, "bulk")[0];
      pred_react = getTemperatures(call_up, id_call, "bulk")[1];

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
std::vector<Real>
ADComputeIntPValRDXMISTERnetNSFull::interpolation(const std::vector<Real> A, const std::vector<Real> B, const Real t){
  std::vector<Real> res;
  res.reserve(A.size());
  for (size_t i = 0; i < A.size(); ++i){
    res.push_back((1. - t) * A[i] + t * B[i]);
  }
  return res;
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

std::vector<Real>
ADComputeIntPValRDXMISTERnetNSFull::getTemperatures(const Real up, const int id, const std::string phase){

  Real lower_bound;
  Real upper_bound;

  Real interval_number;

  std::vector<Real> lower_temps_shock;
  std::vector<Real> upper_temps_shock;

  std::vector<Real> lower_temps_react;
  std::vector<Real> upper_temps_react;

  //here we need to branch based on the phase
  if (phase == "binder"){
    for (unsigned int i = 1; i < _up_values.size(); ++i){
      if (_up_values[i] > up){
        lower_bound = _up_values[i - 1];
        upper_bound = _up_values[i];
        interval_number = i - 1;
        _interval = interval_number;
        
        lower_temps_shock = _temperature_values_shock[i - 1];
        upper_temps_shock = _temperature_values_shock[i];

        lower_temps_react = _temperature_values_react[i - 1];
        upper_temps_react = _temperature_values_react[i];
        break;
      }
    }
  }
  else if (phase == "pore"){
    for (unsigned int j = 1; j < _up_values_pore.size(); ++j){
      if (_up_values_pore[j] > up){
        lower_bound = _up_values_pore[j - 1];
        upper_bound = _up_values_pore[j];
        interval_number = j - 1;
        _interval = interval_number;

        lower_temps_shock = _temperature_values_shock_pore[j - 1];
        upper_temps_shock = _temperature_values_shock_pore[j];

        lower_temps_react = _temperature_values_react_pore[j - 1];
        upper_temps_react = _temperature_values_react_pore[j];
        break;
      }
    }
  }
  else if (phase == "bulk"){
    for (unsigned int k = 1; k < _up_values.size(); ++k){
      if (_up_values[k] > up){
        lower_bound = _up_values[k - 1];
        upper_bound = _up_values[k];
        interval_number = k - 1;
        _interval = interval_number;
        
        lower_temps_shock = _temperature_values_shock[k - 1];
        upper_temps_shock = _temperature_values_shock[k];

        lower_temps_react = _temperature_values_react[k - 1];
        upper_temps_react = _temperature_values_react[k];
        break;
      }
    }
  }

  Real t = (up - lower_bound) / (upper_bound - lower_bound);

  //test: use this same t to interpolate time

  _ratio = t;

  //interpolate temperatures based on t
  //these interpolation steps are naive to the microstructure type
  std::vector<Real> interpolated_temps_shock = interpolation(lower_temps_shock, upper_temps_shock, t); //this had an error
  std::vector<Real> interpolated_temps_react = interpolation(lower_temps_react, upper_temps_react, t);

  //here we need to do branching again based on the type of microstructure
  Real temp_shock;
  Real temp_react;

  if (phase == "binder"){ //traditional branch
    temp_shock = interpolated_temps_shock.at(id);
    temp_react = interpolated_temps_react.at(id);
  }
  if (phase == "pore"){
    //the pore microstructures are 101, 102, 103, so we need to subtract to access the data
    int id_pore = static_cast<int>(std::round(id - _range_pore.front()));

    id_pore = std::max(0, std::min(id_pore, static_cast<int>(interpolated_temps_shock.size()) - 1));
    
    //now access pore data
    temp_shock = interpolated_temps_shock.at(id_pore);
    temp_react = interpolated_temps_react.at(id_pore);
  }
  if (phase == "bulk"){
    //bulk will preserve the nanoPBX data, but will be called all at an insensitive microstructure
    int id_bulk = _bulk_sensitivity; //the most insensitive, this can be tuned until bulk data is available
    
    //now access bulk data
    temp_shock = interpolated_temps_shock.at(id_bulk);
    temp_react = interpolated_temps_react.at(id_bulk);
  }

  return {temp_shock, temp_react};
}

//times will be edited in a similar manner to account for different phases

Real
ADComputeIntPValRDXMISTERnetNSFull::getTimes(const Real up, const int id, const std::string phase){

  //create times lower and upper
  std::vector<Real> lower_times;
  std::vector<Real> upper_times;

  //define lower and upper times
  //here we account for the phase of the element

  if (phase == "binder"){
    //traditional branch
    lower_times = _time_values[_interval];
    upper_times = _time_values[_interval + 1];
  }
  if (phase == "pore"){
    //here we use the pore time values
    lower_times = _time_values_pore[_interval];
    upper_times = _time_values_pore[_interval + 1];
  }
  if (phase == "bulk"){
    //preserve base nanoPBX
    lower_times = _time_values[_interval];
    upper_times = _time_values[_interval + 1];
  }

  //interpolate based on previously computed interval and balance number
  std::vector<Real> interpolated_times = interpolation(lower_times, upper_times, _ratio); //this had an error
  
  //return based on microstructures
  Real time_value;

  if (phase == "binder"){
    time_value = interpolated_times.at(id);
  }

  if (phase == "pore"){
    const int id_pore = id - 101;
    time_value = interpolated_times.at(id_pore);
  }

  if (phase == "bulk"){
    int id_bulk = _bulk_sensitivity;
    time_value = interpolated_times.at(id_bulk);
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