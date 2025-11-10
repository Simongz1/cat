#include "ADComputeIntPValRDXMISTERnetNSFull.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <fstream>

registerMooseObject("mistApp", ADComputeIntPValRDXMISTERnetNSFull);

InputParameters
ADComputeIntPValRDXMISTERnetNSFull::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Standard compute Mie Gruneisen Pressure with JWL pressure for reacted material. Also computes artificial viscosity contribution");
  params.addRequiredParam<Real>("T_ref", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("element_size", "element_size");
  params.addCoupledVar("Y_final", "final products mass fraction");
  params.addRequiredParam<Real>("A_u", "JWL reacted EOS parameter A1");
  params.addRequiredParam<Real>("R1_u", "JWL reacted EOS parameter B1");
  params.addRequiredParam<Real>("B_u", "JWL reacted EOS parameter A2");
  params.addRequiredParam<Real>("R2_u", "JWL reacted EOS parameter B2");
  params.addRequiredParam<Real>("A_r", "JWL reacted EOS parameter A1");
  params.addRequiredParam<Real>("R1_r", "JWL reacted EOS parameter B1");
  params.addRequiredParam<Real>("B_r", "JWL reacted EOS parameter A2");
  params.addRequiredParam<Real>("R2_r", "JWL reacted EOS parameter B2");
  params.addRequiredParam<Real>("omega_u", "JWL reacted EOS parameter omega");
  params.addRequiredParam<Real>("omega_r", "JWL reacted EOS parameter omega");
  params.addRequiredParam<Real>("flag_threshold", "pressure flag threshold value");
  params.addRequiredParam<Real>("P0", "pressure at ambient temperature");
  params.addRequiredParam<Real>("use_RDX", "use RDX model or HMX");
  //request HMX parameters
  params.addRequiredParam<Real>("Gamma", "Grunseisen Parameter");
  params.addRequiredParam<Real>("slope", "slope");
  params.addRequiredParam<Real>("A1", "JWL reacted EOS parameter A1");
  params.addRequiredParam<Real>("R1", "JWL reacted EOS parameter B1");
  params.addRequiredParam<Real>("A2", "JWL reacted EOS parameter A2");
  params.addRequiredParam<Real>("R2", "JWL reacted EOS parameter B2");
  params.addRequiredParam<Real>("omega", "JWL reacted EOS parameter omega");
  //velocity for calling mistnet
  params.addCoupledVar("vx", "x component of velocity");
  params.addCoupledVar("ax", "x component of acceleration");

  //test: using both components of v and a to define shock call
  params.addCoupledVar("vy", "y component of velocity");
  params.addCoupledVar("ay", "y component of acceleration");
  params.addCoupledVar("v_vect", "vector variable that stores velocity components");
  params.addCoupledVar("a_vect", "vector variable that stores acceleration components");
  params.addRequiredParam<Real>("thr_a", "acceleration threshold");
  params.addRequiredParam<Real>("thr_v", "velocity threshold");

  //test: use gradient to compute activation
  

  params.addRequiredParam<Real>("up", "Piston velocity");
  params.addRequiredCoupledVar("density_i", "Coupled value");
  params.addRequiredParam<Real>("mask_size", "mask_size");
  params.addRequiredParam<Real>("use_mask", "use_mask");

  params.addParam<Real>("extrapolate", "extrapolate");
  //CSV
  params.addRequiredParam<std::string>("csv_shock", "csv_shock");
  params.addRequiredParam<std::string>("csv_react", "csv_react");
  params.addRequiredParam<std::string>("csv_times", "csv_times");

  params.addRequiredParam<std::string>("csv_shock_pore", "csv_shock_pore");
  params.addRequiredParam<std::string>("csv_react_pore", "csv_react_pore");
  params.addRequiredParam<std::string>("csv_times_pore", "csv_times_pore");

  params.addRequiredParam<bool>("use_fitted_eos", "use_fitted_eos");
  params.addRequiredParam<bool>("use_EOS_table", "use_EOS_table");
  params.addRequiredParam<bool>("use_magnitude", "use_magnitude");
  params.addRequiredParam<std::string>("csv_unreacted", "csv_unreacted");
  params.addRequiredParam<std::string>("csv_reacted", "csv_reacted");
  params.addRequiredCoupledVar("density_csv", "density_csv");

  //test
  params.addRequiredParam<bool>("use_av_tracking", "use_av_tracking");

  //retrieve bulk ID
  params.addRequiredParam<unsigned int>("bulk_MicroID", "bulk_MicroID");
  params.addParam<unsigned int>("bulk_sensitivity", 20, "bulk_sensitivity");

  //retrieve range for pore IDS
  params.addRequiredParam<std::vector<unsigned int>>("range_pore", "range_pore");
  return params;
}

ADComputeIntPValRDXMISTERnetNSFull::ADComputeIntPValRDXMISTERnetNSFull(
    const InputParameters & parameters)
  : Material(parameters),
  //rest
    _T_ref(getParam<Real>("T_ref")),
    _rho(getADMaterialProperty<Real>("density")),
    _Cv(getADMaterialProperty<Real>("specific_heat")),
    _T(adCoupledValue("temperature")),

    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),

    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _Le(getParam<Real>("element_size")),

    _pressure_mg(declareADProperty<Real>("pressure_mg")),
    _pressure_JWL(declareADProperty<Real>("pressure_JWL")),
    _pressure_total(declareADProperty<Real>("pressure_total")),

    _dP_dT(declareADProperty<Real>("dP_dT")),

    _Y_final(adCoupledValue("Y_final")),

    _A_u(getParam<Real>("A_u")),
    _R1_u(getParam<Real>("R1_u")),
    _B_u(getParam<Real>("B_u")),
    _R2_u(getParam<Real>("R2_u")),
    _omega_u(getParam<Real>("omega_u")),
    _A_r(getParam<Real>("A_r")),
    _R1_r(getParam<Real>("R1_r")),
    _B_r(getParam<Real>("B_r")),
    _R2_r(getParam<Real>("R2_r")),
    _omega_r(getParam<Real>("omega_r")),

    _F(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _F_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),

    _P0(getParam<Real>("P0")),
    _use_RDX(getParam<Real>("use_RDX")),

    //get HMX parameters
    _Gamma(getParam<Real>("Gamma")),
    _s(getParam<Real>("slope")),
    _A1(getParam<Real>("A1")),
    _R1(getParam<Real>("R1")),
    _A2(getParam<Real>("A2")),
    _R2(getParam<Real>("R2")),
    _omega(getParam<Real>("omega")),

    _vx(adCoupledValue("vx")),
    _ax(adCoupledValue("ax")),
    _vy(adCoupledValue("vy")),
    _ay(adCoupledValue("ay")),
    //test:vector variable
    _v_vect(coupledVectorValue("v_vect")),
    _a_vect(coupledVectorValue("a_vect")),

    _thr_a(getParam<Real>("thr_a")),
    _thr_v(getParam<Real>("thr_v")),
    _v_flag(declareProperty<Real>("v_flag")),
    _v_flag_old(getMaterialPropertyOld<Real>("v_flag")),
    _up(getParam<Real>("up")),
    _temperature_mister_shock(declareProperty<Real>("temperature_mister_shock")),
    _temperature_mister_shock_old(getMaterialPropertyOld<Real>("temperature_mister_shock")),
    _temperature_mister_react(declareProperty<Real>("temperature_mister_react")),
    _temperature_mister_react_old(getMaterialPropertyOld<Real>("temperature_mister_react")),
    _density_i(coupledValue("density_i")),
    _called_up(declareProperty<Real>("called_up")),
    _called_up_old(getMaterialPropertyOld<Real>("called_up")),
    _us(declareADProperty<Real>("us")),
    _pressure_av(declareADProperty<Real>("pressure_av")),
    _pressure_av_old(getMaterialPropertyOld<Real>("pressure_av")),

    //get older
    _pressure_av_older(getMaterialPropertyOlder<Real>("pressure_av")),

    _csv_shock(getParam<std::string>("csv_shock")),
    _csv_react(getParam<std::string>("csv_react")),
    _csv_times(getParam<std::string>("csv_times")),

    //csv retrieval for pore
    _csv_shock_pore(getParam<std::string>("csv_shock_pore")),
    _csv_react_pore(getParam<std::string>("csv_react_pore")),
    _csv_times_pore(getParam<std::string>("csv_times_pore")),

    //declare time

    _time_react(declareProperty<Real>("time_react")),
    _time_react_old(getMaterialPropertyOld<Real>("time_react")),
    _use_fitted_eos(getParam<bool>("use_fitted_eos")),
    _use_EOS_table(getParam<bool>("use_EOS_table")),
    _use_magnitude(getParam<bool>("use_magnitude")),
    _csv_unreacted(getParam<std::string>("csv_unreacted")),
    _csv_reacted(getParam<std::string>("csv_reacted")),
    _density(declareADProperty<Real>("density")),
    _density_csv(coupledValue("density_csv")),
    _use_av_tracking(getParam<bool>("use_av_tracking")),

    //bulk grains assignment
    _bulk_MicroID(getParam<unsigned int>("bulk_MicroID")),
    _bulk_sensitivity(getParam<unsigned int>("bulk_sensitivity")),
    _range_pore(getParam<std::vector<unsigned int>>("range_pore"))
{
  //here I cache the table only once
  _csv_total_shock = readCSV(_csv_shock);
  _csv_total_react = readCSV(_csv_react);
  _csv_total_times = readCSV(_csv_times);

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

  for (auto &row : _csv_total_times){
    if (row.size() < 2){
      mooseError("need bigger CSV");
    }
    std::vector<Real> times(row.begin() + 1, row.end());
    _time_values.push_back(std::move(times));
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
  for (auto &row : _csv_total_pu){
    if (row.size() < 2){
      mooseError("need bigger CSV for pressures unreacted");
    }
    _Ju_values.push_back(row[0]);
    _Pu_values.push_back(row[1]);
  }

  for (auto &row : _csv_total_pr){
    if (row.size() < 2){
      mooseError("need bigger CSV for pressures reacted");
    }
    _Jr_values.push_back(row[0]);
    _Pr_values.push_back(row[1]);
  }
  ////////////////////////////////////////////////////////////////
}

void 
ADComputeIntPValRDXMISTERnetNSFull::initQpStatefulProperties()
{
  _v_flag[_qp] = 0.0;
  _stored_shock = _temperature_mister_shock_old[_qp];
  _stored_react = _temperature_mister_react_old[_qp];
  _stored_time  = _time_react_old[_qp];
  _called_up[_qp] = _called_up_old[_qp];
}

void
ADComputeIntPValRDXMISTERnetNSFull::computeQpProperties()
{
  //store density in material property
  _density[_qp] = ADReal(_density_csv[_qp]);
  //KEEP THIS ORDER

  //ORDER 1
  //keep stateful value 1
  _v_flag[_qp] = _v_flag_old[_qp];
  _temperature_mister_shock[_qp] = _temperature_mister_shock_old[_qp];
  _temperature_mister_react[_qp] = _temperature_mister_react_old[_qp];
  _time_react[_qp] = _time_react_old[_qp];
  _called_up[_qp] = _called_up_old[_qp];

  //ORDER 2
  //activate shock heat: call when velocity is bigger than a value 1

  //test: track shock front using artificial viscosity value

  Real condition_v;
  Real condition_a;
  
  if(_use_magnitude){
    condition_v = _v_vect[_qp].norm();
    condition_a = _a_vect[_qp].norm();

  }else{
    condition_v = std::abs(MetaPhysicL::raw_value(_vx[_qp]));
    condition_a = std::abs(MetaPhysicL::raw_value(_ax[_qp]));
  }

  //////compute AV here
  const RankTwoTensor I2(RankTwoTensor::initIdentity);
  ADReal P_av;
  const Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  ADReal ss = std::sqrt(K0 / _rho[_qp]);

  Real Je = _F[_qp].det();
  Real Je_dot = (_F[_qp].det() - _F_old[_qp].det()) / _dt;

  P_av = _C0 * _rho[_qp] * (Je_dot * std::abs(Je_dot) / std::pow(Je, 2.0)) * std::pow(_Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (Je_dot / Je) * _Le;

  bool _call_condition = false;

  if (_use_av_tracking){
    Real av = MetaPhysicL::raw_value(P_av);
    Real av_old = std::max(_pressure_av_old[_qp], 0.);
    Real av_older = std::max(_pressure_av_older[_qp], 0.);

    //previous timestep derivative
    Real av_diff = (av - av_old) / _dt;
    Real av_diff_old = (av_old - av_older) / _dt;
    Real av_2diff = (av_diff - av_diff_old) / (2 * _dt);
    
    if (av_older < 0. && std::abs(av_old) < 1e-1 && av_diff > 0.){
      _call_condition = true;
    }
  }
  else{
    _call_condition = (_qp == 0. && _v_flag[_qp] == 0. && condition_v > _thr_v && condition_a < _thr_a) ? true : false;
  }

  if (_call_condition){ //we need V and A constraints to make sure the call happens at the actual shock velocity
    _v_flag[_qp] = 1.0; //set flag to one, call in misternet material as condition for ComputeQpProperties
    //get temperatures
    const int id_call = static_cast<int>(std::round(_density_i[_qp]));

    //test: clamp velocity to values
    //here we need to branch if we are using bulk or MicroID values

    //declare values for binder
    Real pred_shock, pred_react, pred_time;

    //define cases
  
    const bool is_pore = (id_call >= _range_pore[0] && id_call <= _range_pore[1]);
    const bool is_bulk = (id_call == _bulk_MicroID);
    const bool is_binder = (!is_pore && !is_bulk);

    //define call velocity explicitly
    const Real call_up = std::clamp(std::abs(MetaPhysicL::raw_value(_vx[_qp])), 0.0, 4.89);

    //retrieve data based on grain or binder
    if (is_binder){ //this is the usual loop
      //assign hand coded values as placeholders
      pred_shock = getTemperatures(call_up, id_call, "binder")[0];
      pred_react = getTemperatures(call_up, id_call, "binder")[1];
      pred_time  = getTimes(call_up, id_call, "binder");
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
      pred_time  = getTimes(call_up, id_call, "bulk");
    }
  
    //store in material property
    _temperature_mister_shock[_qp] = pred_shock;
    _temperature_mister_react[_qp] = pred_react;
    _time_react[_qp]               = pred_time;

    //get the temeprature at the initial qp
    _stored_shock = pred_shock;
    _stored_react = pred_react;
    _stored_time  = pred_time;

    //store the called up value
    _called_up[_qp] = std::abs(MetaPhysicL::raw_value(_vx[_qp]));
  }

  //ORDER 3

  if (_qp > 0.){
    _v_flag[_qp] = _v_flag[0];
    _temperature_mister_shock[_qp] = _temperature_mister_shock[0];
    _temperature_mister_react[_qp] = _temperature_mister_react[0];
    _time_react[_qp]               = _time_react[0];
    _called_up[_qp] = _called_up[0];
  }

  //initialize the symmetric identity tensors

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour
  
  //Compute Mie Gruneisen pressure for unreacted material
  ADReal P_mg;
  ADReal P_JWL;
  //this should only depend on the elastic deformation gradient
  
  Real eta = 1. - Je;
  
  P_mg = _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / Je); //initial thermal expansion term
  P_mg += K0 * eta * (1.0 - (_Gamma / 2.0) * (eta)) / std::pow((1.0 - _s * eta), 2.0);
  
  if(_use_fitted_eos){
    P_mg = _A_u * std::exp(- _R1_u * Je) + _B_u * std::exp(- _R2_u * Je);
    //test: don't use reference temperature
    P_mg += _omega_u * _rho[_qp] * _Cv[_qp] * (_T[_qp]) / Je;
  }
  _pressure_mg[_qp] = - P_mg; //store pressure 

  P_JWL = _A1 * (1.0 - _omega / (_R1 * Je)) * std::exp(- _R1 * Je); //mechanical term 1
  P_JWL += _A2 * (1.0 - _omega / (_R2 * Je)) * std::exp(- _R2 * Je); //mechanical term 2
  P_JWL += _omega * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) / Je; //thermal expansion term

  if(_use_fitted_eos){
    P_JWL = _A_r * std::exp(- _R1_r * Je) + _B_r * std::exp(- _R2_r * Je);
    //test: don't use reference temperature
    P_JWL += _omega_r * _rho[_qp] * _Cv[_qp] * (_T[_qp]) / Je;
  }
  _pressure_JWL[_qp] = - P_JWL; //store

  //pressure interpolation
  ADReal P_total;
  P_total = ((1.0 - _Y_final[_qp]) * (P_mg)) + (_Y_final[_qp] * P_JWL); //compute total pressure
  _pressure_total[_qp] = - P_total;

  //define derivatives of pressure with respect to temperature
  ADReal dPmg_dT;
  ADReal dPJWL_dT;

  dPmg_dT = _Gamma * _rho[_qp] * _Cv[_qp] * (1. / Je);
  dPJWL_dT = _omega * _rho[_qp] * _Cv[_qp] * (1. / Je);
  
  if (_use_fitted_eos){
    dPmg_dT = _omega_u * _rho[_qp] * _Cv[_qp] * (1. / Je);
    dPJWL_dT = _omega_r * _rho[_qp] * _Cv[_qp] * (1. / Je);
  }

  //test: use table

  if(_use_EOS_table){ 
    Real Jac = Je; //elastic volume change, should accont for thermal expansion too
    std::vector<Real> pressures = getPressures(Je);
    _pressure_mg[_qp] = - pressures[0];
    _pressure_JWL[_qp] = - pressures[1];
    _pressure_total[_qp] = - ((1. - _Y_final[_qp]) * pressures[0] + _Y_final[_qp] * pressures[1]);
  }

  //test: use AD to get the derivative of pressure WRT temperature, then write into an ADMaterialProperty
  _dP_dT[_qp] = (1. - _Y_final[_qp]) * dPmg_dT + _Y_final[_qp] * dPJWL_dT;

  //include artificial viscosity
  //ADReal P_av;
  //P_av = _C0 * _rho[_qp] * (Je_dot * std::abs(Je_dot) / std::pow(Je, 2.0)) * std::pow(_Le, 2.0);
  //P_av += _C1 * _rho[_qp] * ss * (Je_dot / Je) * _Le;
  _pressure_av[_qp] = P_av;

  _pressure_total[_qp] +=  P_av;

  //write total pressure into stress as a hydrostatic component
  //compute and declare the derivatives of each partial pressure wrt temperature to consume on PressureHS

  //compute shock velocity
  ADReal us;
  
  if (_use_fitted_eos){
    us = 4.0790 + 1.9370 * _v_vect[_qp].norm();
  }else{
    us = ss + (_s * _v_vect[_qp].norm());
  }
  _us[_qp] = us;
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
    const int id_bulk = 0;
    time_value = interpolated_times.at(id_bulk);
  }

  //return whatever value went through
  return time_value;
}

//std::vector<Real>
//ADComputeIntPValRDXMISTERnetNSFull::getPressures(const Real Jac){
//  //this gets the unreacted pressure for a given J = det(F) from the csv table
//  Real P;
//  Real Pu;
//  Real Pr;
//  Real lower_Ju;
//  Real upper_Ju;
//  Real lower_Jr;
//  Real upper_Jr;
//  Real lower_Pu;
//  Real upper_Pu;
//  Real lower_Pr;
//  Real upper_Pr;
//
//  //run thorugh J values to get the interval where the actual J is
//  Real loc_u;
//  Real clamped_Ju = std::clamp(Jac, _Ju_values[0], _Ju_values.back());
//  for (unsigned int i = 1; i < _Ju_values.size(); ++i){
//    if(_Ju_values[i] > clamped_Ju){
//      lower_Ju = _Ju_values[i - 1];
//      upper_Ju = _Ju_values[i];
//      loc_u = i;
//      break;
//    }
//  }
//  Real tu = (clamped_Ju - lower_Ju) / (upper_Ju - lower_Ju);
//
//  Real loc_r;
//  Real clamped_Jr = std::clamp(Jac, _Jr_values[0], _Jr_values.back());
//  for (unsigned int i = 1; i < _Jr_values.size(); ++i){
//    if(_Jr_values[i] > clamped_Jr){ //make sure we don't go out of bounds, specially since data for reacted is limited
//      lower_Jr = _Jr_values[i - 1];
//      upper_Jr = _Jr_values[i];
//      loc_r = i;
//      break;
//    }
//  }
//  Real tr = (clamped_Jr - lower_Jr) / (upper_Jr - lower_Jr);
//  
//  //use this interval to intepolate the pressures
//  lower_Pu = _Pu_values[loc_u - 1];
//  upper_Pu = _Pu_values[loc_u];
//  lower_Pr = _Pr_values[loc_r - 1];
//  upper_Pr = _Pr_values[loc_r];
//
//  Pu = (1. - tu) * lower_Pu + tu * upper_Pu;
//  Pr = (1. - tr) * lower_Pr + tr * upper_Pr;
//  return {Pu, Pr};
//}

std::vector<Real>
ADComputeIntPValRDXMISTERnetNSFull::getPressures(const Real J)
{
  auto interp_no_extrap = [](Real x,
                             const std::vector<Real> & X,
                             const std::vector<Real> & Y) -> Real
  {
    auto it = std::lower_bound(X.begin(), X.end(), x);

    if (it == X.begin())
      return Y.front();
    if (it == X.end())
      return Y.back();

    const size_t i = static_cast<size_t>(std::distance(X.begin(), it));
    const Real x0 = X[i - 1], x1 = X[i];
    const Real y0 = Y[i - 1], y1 = Y[i];

    const Real dx = x1 - x0;
    if (dx == 0.0)
      return y1;

    const Real t = (x - x0) / dx;
    
    return (1.0 - t) * y0 + t * y1;
  };

  const Real Pu = interp_no_extrap(J, _Ju_values, _Pu_values);
  const Real Pr = interp_no_extrap(J, _Jr_values, _Pr_values);
  return {Pu, Pr};
}

//test: generate a function to retrieve temperatures, times, and densities
