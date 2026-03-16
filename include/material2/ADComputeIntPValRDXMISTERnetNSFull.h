
#pragma once

#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "MathUtils.h"
#include "RankFourTensor.h"
#include <vector>
#include "ElasticityTensorTools.h"
#include "MooseUtils.h"
#include "DistributionInterface.h"

class ADComputeIntPValRDXMISTERnetNSFull : public Material
{
public:
  static InputParameters validParams();

  ADComputeIntPValRDXMISTERnetNSFull(const InputParameters & parameters);

//  void initialSetup() override;

protected:


  const std::vector<VariableName> _v_components;
  const std::vector<VariableName> _a_components;

  const Real _thr_a;
  const Real _thr_v;
  ADMaterialProperty<Real> &_v_flag;
  const MaterialProperty<Real> &_v_flag_old;

  MaterialProperty<Real> & _temperature_mister_shock;
  const MaterialProperty<Real> &_temperature_mister_shock_old;
  MaterialProperty<Real> & _temperature_mister_react;
  const MaterialProperty<Real> &_temperature_mister_react_old;

  const VariableValue & _density_i;

  Real _stored_shock;
  Real _stored_react;
  Real _stored_time;

  ADMaterialProperty<Real> &_called_up;
  const MaterialProperty<Real> &_called_up_old;
  ADMaterialProperty<Real> &_us;

  bool _call_condition;

  // csv
  const std::string _csv_shock;
  const std::string _csv_react;
  const std::string _csv_times;

  const std::string _csv_shock_pore;
  const std::string _csv_react_pore;
  const std::string _csv_times_pore;

  ADMaterialProperty<Real> &_time_react;
  const MaterialProperty<Real> &_time_react_old;

  const std::string _csv_unreacted;
  const std::string _csv_reacted;
  ADMaterialProperty<Real> &_density;

  const VariableValue &_density_csv;
  const unsigned int _bulk_MicroID;
  const unsigned int _bulk_sensitivity;
  const std::vector<unsigned int> _range_pore;

  ADMaterialProperty<Real> &_rate_tracking;
  const Real _h;
  const ADVariableValue &_tracking;
  const bool _use_tabular_time;
  const bool _use_distributions;
  
  Distribution const *_distribution_lower;
  Distribution const *_distribution_upper;

  /////////////////////////////////

  std::vector<const ADVariableValue *> _v;
  std::vector<const ADVariableValue *> _a;

  std::vector<std::vector<Real>> _csv_total_shock;
  std::vector<std::vector<Real>> _csv_total_react;
  std::vector<std::vector<Real>> _csv_total_times;


  std::vector<std::vector<Real>> _csv_total_shock_pore;
  std::vector<std::vector<Real>> _csv_total_react_pore;
  std::vector<std::vector<Real>> _csv_total_times_pore;

  //
  std::vector<std::vector<Real>> _csv_total_pu;
  std::vector<std::vector<Real>> _csv_total_pr;
  std::vector<Real> _up_values;
  std::vector<Real> _up_values_pore;
  std::vector<std::vector<Real>> _temperature_values_shock;
  std::vector<std::vector<Real>> _temperature_values_react;
  std::vector<std::vector<Real>> _temperature_values_shock_pore;
  std::vector<std::vector<Real>> _temperature_values_react_pore;
  std::vector<std::vector<Real>> _time_values;
  std::vector<std::vector<Real>> _time_values_pore;

  //definitions for global interpolation

  Real _interval;
  Real _ratio;

  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  virtual std::vector<Real> interpolation(const std::vector<Real> A, const std::vector<Real> B, const Real t);
  virtual std::vector<std::vector<Real>> readCSV(const std::string csv_name);
  virtual std::vector<Real> getTemperatures(const Real up, const int id, const std::string phase);
  virtual Real getTimes(const Real up, const int id, const std::string phase);
  virtual Real getDistributionTime(const Real up, const Real predicted_temp);
};
