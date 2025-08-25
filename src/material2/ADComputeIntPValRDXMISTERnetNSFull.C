#include "ADComputeIntPValRDXMISTERnetNSFull.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <fstream>

registerMooseObject("TensorMechanicsApp", ADComputeIntPValRDXMISTERnetNSFull);

InputParameters
ADComputeIntPValRDXMISTERnetNSFull::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Standard compute Mie Gruneisen Pressure with JWL pressure for reacted material. Also computes artificial viscosity contribution");
  params.addRequiredParam<Real>("T_ref", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("C0", "artificial viscosity C0 parameter");
  params.addRequiredParam<Real>("C1", "artificial viscosity C1 parameter");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredParam<Real>("element_size", "element_size");
  params.addRequiredCoupledVar("Y_final", "final products mass fraction");
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
  params.addRequiredCoupledVar("vx", "x component of velocity");
  params.addRequiredCoupledVar("ax", "x component of acceleration");
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

  params.addRequiredParam<bool>("use_fitted_eos", "use_fitted_eos");
  return params;
}

ADComputeIntPValRDXMISTERnetNSFull::ADComputeIntPValRDXMISTERnetNSFull(
    const InputParameters & parameters)
  : Material(parameters),
  //rest
    _T_ref(getParam<Real>("T_ref")),
    _rho(getADMaterialProperty<Real>("density")),
    _Cv(getADMaterialProperty<Real>("specific_heat")),
    _T(coupledValue("temperature")),

    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),

    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _Le(getParam<Real>("element_size")),

    _pressure_mg(declareADProperty<Real>("pressure_mg")),
    _pressure_JWL(declareADProperty<Real>("pressure_JWL")),
    _pressure_total(declareADProperty<Real>("pressure_total")),

    _dP_dT(declareADProperty<Real>("dP_dT")),

    _Y_final(coupledValue("Y_final")),

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

    _vx(coupledValue("vx")),
    _vx_old(coupledValueOld("vx")),
    _ax(coupledValue("ax")),

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
    _mask_size(getParam<Real>("mask_size")),
    _use_mask(getParam<Real>("use_mask")),
    _called_up(declareProperty<Real>("called_up")),
    _called_up_old(getMaterialPropertyOld<Real>("called_up")),
    _us(declareADProperty<Real>("us")),
    _pressure_av(declareADProperty<Real>("pressure_av")),

    _csv_shock(getParam<std::string>("csv_shock")),
    _csv_react(getParam<std::string>("csv_react")),
    _csv_times(getParam<std::string>("csv_times")),

    //declare time

    _time_react(declareProperty<Real>("time_react")),
    _time_react_old(getMaterialPropertyOld<Real>("time_react")),
    _use_fitted_eos(getParam<bool>("use_fitted_eos"))
  //here I cache the table only once

{
  //here I cache the table only once
  _csv_total_shock = readCSV(_csv_shock);
  _csv_total_react = readCSV(_csv_react);
  _csv_total_times = readCSV(_csv_times);

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

  //test: get times into usable array

  for (auto &row : _csv_total_times){
    if (row.size() < 2){
      mooseError("need bigger CSV");
    }
    std::vector<Real> times(row.begin() + 1, row.end());
    _time_values.push_back(std::move(times));
  }
}

void 
ADComputeIntPValRDXMISTERnetNSFull::initQpStatefulProperties()
{
  _v_flag[_qp] = 0.0; //initialize flag
  _stored_shock = _temperature_mister_shock_old[_qp];
  _stored_react = _temperature_mister_react_old[_qp];
  _stored_time  = _time_react_old[_qp];
  _called_up[_qp] = _called_up_old[_qp];

  ///////////////////////////
}

void
ADComputeIntPValRDXMISTERnetNSFull::computeQpProperties()
{

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

  if (_qp == 0. && _v_flag[_qp] == 0. && std::abs(_vx[_qp]) > _thr_v && std::abs(_ax[_qp]) < _thr_a){ //we need V and A constraints to make sure the call happens at the actual shock velocity
    _v_flag[_qp] = 1.0; //set flag to one, call in misternet material as condition for ComputeQpProperties
    //get temperatures
    const int id_call = (_density_i[_qp]);

    Real pred_shock = getTemperatures(std::abs(_vx[_qp]), id_call)[0];
    Real pred_react = getTemperatures(std::abs(_vx[_qp]), id_call)[1];
    Real pred_time  = getTimes(std::abs(_vx[_qp]), id_call);

    //store in material property
    _temperature_mister_shock[_qp] = pred_shock;
    _temperature_mister_react[_qp] = pred_react;
    _time_react[_qp]               = pred_time;

    //get the temeprature at the initial qp
    _stored_shock = pred_shock;
    _stored_react = pred_react;
    _stored_time  = pred_time;

    //store the called up value
    _called_up[_qp] = std::abs(_vx[_qp]);
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
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour
  Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  ADReal ss = std::sqrt(K0 / _rho[_qp]);

  //Compute Mie Gruneisen pressure for unreacted material
  ADReal P_mg;
  ADReal P_JWL;
  //this should only depend on the elastic deformation gradient
  Real Je = _F[_qp].det();
  Real Je_dot = (_F[_qp].det() - _F_old[_qp].det()) / _dt;
  Real eta = 1. - Je;
  
  P_mg = _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / Je); //initial thermal expansion term
  P_mg += K0 * eta * (1.0 - (_Gamma / 2.0) * (eta)) / std::pow((1.0 - _s * eta), 2.0);
  
  if(_use_fitted_eos){
    P_mg = _A_u * std::exp(- _R1_u * Je) + _B_u * std::exp(- _R2_u * Je);
    P_mg += _omega_u * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) / Je;
  }
  _pressure_mg[_qp] = - P_mg; //store pressure 

  P_JWL = _A1 * (1.0 - _omega / (_R1 * Je)) * std::exp(- _R1 * Je); //mechanical term 1
  P_JWL += _A2 * (1.0 - _omega / (_R2 * Je)) * std::exp(- _R2 * Je); //mechanical term 2
  P_JWL += _omega * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) / Je; //thermal expansion term

  if(_use_fitted_eos){
    P_JWL = _A_r * std::exp(- _R1_r * Je) + _B_r * std::exp(- _R2_r * Je);
    P_JWL += _omega_r * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) / Je;
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

  _dP_dT[_qp] = ((1.0 - _Y_final[_qp]) * (dPmg_dT)) + (_Y_final[_qp] * dPJWL_dT);

  //include artificial viscosity
  ADReal P_av;
  P_av = _C0 * _rho[_qp] * (Je_dot * std::abs(Je_dot) / std::pow(Je, 2.0)) * std::pow(_Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (Je_dot / Je) * _Le;
  _pressure_av[_qp] = P_av;

  _pressure_total[_qp] +=  P_av;

  //write total pressure into stress as a hydrostatic component
  //compute and declare the derivatives of each partial pressure wrt temperature to consume on PressureHS

  //compute shock velocity
  _us[_qp] = ss + (_s * std::abs(_vx[_qp]));
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
ADComputeIntPValRDXMISTERnetNSFull::getTemperatures(const Real up, const int id){

  Real lower_bound;
  Real upper_bound;
  Real interval_number;
  std::vector<Real> lower_temps_shock;
  std::vector<Real> upper_temps_shock;

  std::vector<Real> lower_temps_react;
  std::vector<Real> upper_temps_react;

  for (unsigned int i = 1; i < _up_values.size(); ++i){
    if (_up_values[i] > up) {
      lower_bound = _up_values[i - 1];
      upper_bound = _up_values[i];
      interval_number = i - 1;

      //test: save the interval number on the fly
      _interval = interval_number;

      lower_temps_shock = _temperature_values_shock[i - 1];
      upper_temps_shock = _temperature_values_shock[i];

      lower_temps_react = _temperature_values_react[i - 1];
      upper_temps_react = _temperature_values_react[i];
      break;
    }
  }
  

  Real t = (up - lower_bound) / (upper_bound - lower_bound);

  //test: use this same t to interpolate time

  _ratio = t;

  //interpolate temperatures based on t

  std::vector<Real> interpolated_temps_shock = interpolation(lower_temps_shock, upper_temps_shock, t); //this had an error
  std::vector<Real> interpolated_temps_react = interpolation(lower_temps_react, upper_temps_react, t);
  Real temp_shock = interpolated_temps_shock.at(id);
  Real temp_react = interpolated_temps_react.at(id);

  //add the interpolation of the time to deflagration with the new data
  
  return {temp_shock, temp_react};
}

Real
ADComputeIntPValRDXMISTERnetNSFull::getTimes(const Real up, const int id){

  //create times lower and upper
  std::vector<Real> lower_times;
  std::vector<Real> upper_times;

  //define lower and upper times

  lower_times = _time_values[_interval];
  upper_times = _time_values[_interval + 1];

  //interpolate based on previously computed interval and balance number
  std::vector<Real> interpolated_times = interpolation(lower_times, upper_times, _ratio); //this had an error
  
  return interpolated_times.at(id);
}