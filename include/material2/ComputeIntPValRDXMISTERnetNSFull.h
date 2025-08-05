
#pragma once

#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "MathUtils.h"
#include "RankFourTensor.h"
#include <vector>

class ComputeIntPValRDXMISTERnetNSFull : public Material
{
public:
  static InputParameters validParams();

  ComputeIntPValRDXMISTERnetNSFull(const InputParameters & parameters);

//  void initialSetup() override;

protected:
  const Real _T_ref;
  const MaterialProperty<Real> &_rho;
  const MaterialProperty<Real> &_Cv;
  const VariableValue &_T;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain_old;
  const Real _C0;
  const Real _C1;

  //variable element length computation
  const MaterialProperty<RankFourTensor> &_elasticity_tensor;
  const Real _Le;
  MaterialProperty<Real> &_pressure_mg;
  MaterialProperty<Real> &_pressure_JWL;
  MaterialProperty<Real> &_pressure_total;
  const VariableValue &_Y_final;
  const Real _A_u;
  const Real _R1_u;
  const Real _B_u;
  const Real _R2_u;
  const Real _omega_u;
  const Real _A_r;
  const Real _R1_r;
  const Real _B_r;
  const Real _R2_r;
  const Real _omega_r;
  MaterialProperty<RankTwoTensor> &_dMG_dT;
  MaterialProperty<RankTwoTensor> &_dJWL_dT;
  const MaterialProperty<Real> &_pressure_total_old;
  MaterialProperty<Real> &_pressure_hist;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> &_f_inv;
  const Real _P0;
  const Real _use_RDX;
  const Real _Gamma;
  const Real _s;
  const Real _A1;
  const Real _R1;
  const Real _A2;
  const Real _R2;
  const Real _omega;
  const VariableValue &_vx;
  const VariableValue &_vx_old;
  const VariableValue &_ax;
  const Real _thr_a;
  const Real _thr_v;
  MaterialProperty<Real> &_v_flag;
  const MaterialProperty<Real> &_v_flag_old;

  const Real _up;
  MaterialProperty<Real> & _temperature_mister_shock;
  const MaterialProperty<Real> &_temperature_mister_shock_old;
  MaterialProperty<Real> & _temperature_mister_react;
  const MaterialProperty<Real> &_temperature_mister_react_old;

  const VariableValue & _density_i;
  const Real _mask_size;
  const Real _use_mask;

  Real _stored_shock;
  Real _stored_react;

  MaterialProperty<Real> &_called_up;
  const MaterialProperty<Real> &_called_up_old;
  MaterialProperty<Real> &_us;
  MaterialProperty<RankTwoTensor> &_pressure_av;

  //testing ignition interpolation
  
  const VariableValue &_lambda;
  const Real _use_IG;

  const std::string _csv_shock;
  const std::string _csv_react;

  std::vector<std::vector<Real>> _csv_total_shock;
  std::vector<std::vector<Real>> _csv_total_react;
  std::vector<Real> _up_values;
  std::vector<std::vector<Real>> _temperature_values_shock;
  std::vector<std::vector<Real>> _temperature_values_react;

  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  virtual std::vector<Real> interpolation(const std::vector<Real> A, const std::vector<Real> B, const Real t);
  virtual std::vector<std::vector<Real>> readCSV(const std::string csv_name);
  virtual std::vector<Real> getTemperatures(const Real up, const int id);
};
