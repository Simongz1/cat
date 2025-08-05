#include "ComputeIntPValRDXMISTERnetNSFull.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <fstream>

registerMooseObject("TensorMechanicsApp", ComputeIntPValRDXMISTERnetNSFull);

InputParameters
ComputeIntPValRDXMISTERnetNSFull::validParams()
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

  params.addRequiredParam<Real>("up", "Piston velocity");
  params.addRequiredCoupledVar("density_i", "Coupled value");
  params.addRequiredParam<Real>("mask_size", "mask_size");
  params.addRequiredParam<Real>("use_mask", "use_mask");
  params.addRequiredCoupledVar("lambda", "lambda");
  params.addRequiredParam<Real>("use_IG", "use_IG");

  params.addParam<Real>("extrapolate", "extrapolate");
  //CSV
  params.addRequiredParam<std::string>("csv_shock", "csv_shock");
  params.addRequiredParam<std::string>("csv_react", "csv_react");
  return params;
}

ComputeIntPValRDXMISTERnetNSFull::ComputeIntPValRDXMISTERnetNSFull(
    const InputParameters & parameters)
  : Material(parameters),
  //rest
    _T_ref(getParam<Real>("T_ref")),
    _rho(getMaterialProperty<Real>("density")),
    _Cv(getMaterialProperty<Real>("specific_heat")),
    _T(coupledValue("temperature")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>("mechanical_strain")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
    _Le(getParam<Real>("element_size")),
    _pressure_mg(declareProperty<Real>("pressure_mg")),
    _pressure_JWL(declareProperty<Real>("pressure_JWL")),
    _pressure_total(declareProperty<Real>("pressure_total")),
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
    _dMG_dT(declareProperty<RankTwoTensor>("dMG_dT")),
    _dJWL_dT(declareProperty<RankTwoTensor>("dJWL_dT")),
    _pressure_total_old(getMaterialPropertyOld<Real>("pressure_total")),
    _pressure_hist(declareProperty<Real>("pressure_hist")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _f_inv(getMaterialProperty<RankTwoTensor>("inverse_incremental_deformation_gradient")),
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
    _us(declareProperty<Real>("us")),
    _pressure_av(declareProperty<RankTwoTensor>("pressure_av")),
    //testing ignition interpolation
    _lambda(coupledValue("lambda")),
    _use_IG(getParam<Real>("use_IG")),
    _csv_shock(getParam<std::string>("csv_shock")),
    _csv_react(getParam<std::string>("csv_react"))
{
  //here I cache the table only once
  _csv_total_shock = readCSV(_csv_shock);
  _csv_total_react = readCSV(_csv_react);

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
}

void 
ComputeIntPValRDXMISTERnetNSFull::initQpStatefulProperties()
{
  _v_flag[_qp] = 0.0; //initialize flag
  _stored_shock = _temperature_mister_shock_old[_qp];
  _stored_react = _temperature_mister_react_old[_qp];
  _called_up[_qp] = _called_up_old[_qp];

  ///////////////////////////
}

void
ComputeIntPValRDXMISTERnetNSFull::computeQpProperties()
{

  //KEEP THIS ORDER

  //ORDER 1
  //keep stateful value 1
  _v_flag[_qp] = _v_flag_old[_qp];
  _temperature_mister_shock[_qp] = _temperature_mister_shock_old[_qp];
  _temperature_mister_react[_qp] = _temperature_mister_react_old[_qp];
  _called_up[_qp] = _called_up_old[_qp];

  //ORDER 2
  //activate shock heat: call when velocity is bigger than a value 1
  Real variation;
  variation = std::abs((_vx[_qp] - _vx_old[_qp]) / _dt);
  if (_qp == 0. && _v_flag[_qp] == 0. && std::abs(_vx[_qp]) > _thr_v && std::abs(_ax[_qp]) < _thr_a){ //we need V and A constraints to make sure the call happens at the actual shock velocity
    _v_flag[_qp] = 1.0; //set flag to one, call in misternet material as condition for ComputeQpProperties
    //get temperatures
    const int id_call = (_density_i[_qp]);

    Real pred_shock = getTemperatures(std::abs(_vx[_qp]), id_call)[0];
    Real pred_react = getTemperatures(std::abs(_vx[_qp]), id_call)[1];

    //store in material property
    _temperature_mister_shock[_qp] = pred_shock;
    _temperature_mister_react[_qp] = pred_react;

    //get the temeprature at the initial qp
    _stored_shock = pred_shock;
    _stored_react = pred_react;

    //store the called up value
    _called_up[_qp] = std::abs(_vx[_qp]);
  }

  //ORDER 3

  if (_qp > 0.){
    _v_flag[_qp] = _v_flag[0];
    _temperature_mister_shock[_qp] = _temperature_mister_shock[0];
    _temperature_mister_react[_qp] = _temperature_mister_react[0];
    _called_up[_qp] = _called_up[0];
  }

  //initialize the symmetric identity tensors
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  //compute sound speed and bulk modulus from elasticity tensors
  //this is important for the case later on when we add anisotropic behaviour
  Real K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  Real ss = std::sqrt(K0 / _rho[_qp]);

  //Compute Mie Gruneisen pressure for unreacted material
  Real P_mg;
  Real P_JWL;
  //
  Real J = _deformation_gradient[_qp].det();
  Real J_dot = (_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det()) / _dt;
  Real eta = 1.0 - J;
  Real rho0_rho = 1.0 - eta; //express rho / rho_0
  
  RankTwoTensor epsilon_dot = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;

  P_mg = _Gamma * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) * (1.0 / J); //initial thermal expansion term
  P_mg += K0 * eta * (1.0 - (_Gamma / 2.0) * (eta)) / std::pow((1.0 - _s * eta), 2.0);
  _pressure_mg[_qp] = - P_mg; //store pressure 

  P_JWL = _A1 * (1.0 - _omega / (_R1 * J)) * std::exp(- _R1 * J); //mechanical term 1
  P_JWL += _A2 * (1.0 - _omega / (_R2 * J)) * std::exp(- _R2 * J); //mechanical term 2
  P_JWL += _omega * _rho[_qp] * _Cv[_qp] * (_T[_qp] - _T_ref) / J; //thermal expansion term
  _pressure_JWL[_qp] = - P_JWL; //store

  //pressure interpolation
  Real P_total;

  if (_use_IG == 1.){
    P_total = ((1.0 - _lambda[_qp]) * (P_mg)) + (_lambda[_qp] * P_JWL);
  }

  else {
    P_total = ((1.0 - _Y_final[_qp]) * (P_mg)) + (_Y_final[_qp] * P_JWL); //compute total pressure
  }
        
  _pressure_total[_qp] = - P_total;

  //include artificial viscosity
  Real P_av;
  Real Je;
  Real Je_dot;
  Je = _deformation_gradient[_qp].det();
  Je_dot = ((_deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det()) / _dt);

  P_av = _C0 * _rho[_qp] * (Je_dot * std::abs(Je_dot) / std::pow(Je, 2.0)) * std::pow(_Le, 2.0);
  P_av += _C1 * _rho[_qp] * ss * (Je_dot / Je) * _Le;
  _pressure_av[_qp] = P_av * I2;

  _pressure_total[_qp] += (1. / 3.) * _pressure_av[_qp].trace();

  //write total pressure into stress as a hydrostatic component
  //compute and declare the derivatives of each partial pressure wrt temperature to consume on PressureHS

  //compression heat source declaration
  RankTwoTensor dMG_dT;
  RankTwoTensor dJWL_dT;

  dMG_dT = - _omega_u * _rho[_qp] * _Cv[_qp] * (1.0 / J);
  dJWL_dT = - _omega_r * _rho[_qp] * _Cv[_qp] * (1.0 / J);

  _dMG_dT[_qp] = dMG_dT;
  _dJWL_dT[_qp] = dJWL_dT;

  //compute lagrangian strain and other tensors
  RankTwoTensor f = _f_inv[_qp].inverse();
  RankTwoTensor f_bar = f / std::cbrt(f.det());

  //compute shock velocity
  _us[_qp] = ss + (_s * std::abs(_vx[_qp]));
}

//interpolate between values
std::vector<Real>
ComputeIntPValRDXMISTERnetNSFull::interpolation(const std::vector<Real> A, const std::vector<Real> B, const Real t){
  std::vector<Real> res;
  res.reserve(A.size());
  for (size_t i = 0; i < A.size(); ++i){
    res.push_back((1. - t) * A[i] + t * B[i]);
  }
  return res;
}

//helper function to read CSV file
std::vector<std::vector<Real>>
ComputeIntPValRDXMISTERnetNSFull::readCSV(const std::string csv_file_name){
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
ComputeIntPValRDXMISTERnetNSFull::getTemperatures(const Real up, const int id){

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
      lower_temps_shock = _temperature_values_shock[i - 1];
      upper_temps_shock = _temperature_values_shock[i];

      lower_temps_react = _temperature_values_react[i - 1];
      upper_temps_react = _temperature_values_react[i];
      break;
    }
  }
  

  Real t = (up - lower_bound) / (upper_bound - lower_bound);

  //interpolate temperatures based on t

  std::vector<Real> interpolated_temps_shock = interpolation(lower_temps_shock, upper_temps_shock, t); //this had an error
  std::vector<Real> interpolated_temps_react = interpolation(lower_temps_react, upper_temps_react, t);
  Real temp_shock = interpolated_temps_shock.at(id);
  Real temp_react = interpolated_temps_react.at(id);
  
  return {temp_shock, temp_react};
}