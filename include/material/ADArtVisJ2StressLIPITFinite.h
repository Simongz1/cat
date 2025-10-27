#pragma once

#include "ComputeLagrangianStressPK1.h"
#include "GuaranteeConsumer.h"
#include "ElasticityTensorTools.h"
#include "SingleVariableReturnMappingSolution.h"
#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADArtVisJ2StressLIPITFinite
  : public DerivativeMaterialInterface<ComputeLagrangianStressPK1>,
    public GuaranteeConsumer,
    public SingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADArtVisJ2StressLIPITFinite(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeQpPK1Stress() override;

  /// @{ The return mapping residual and derivative
  virtual Real computeReferenceResidual(const Real & effective_trial_stress,
                                        const Real & scalar) override;
  virtual Real computeResidual(const Real & effective_trial_stress, const Real & scalar) override;
  virtual Real computeDerivative(const Real & effective_trial_stress, const Real & scalar) override;
  virtual void
  preStep(const Real & scalar_old, const Real & residual, const Real & jacobian) override;
  /// @}

  const MaterialPropertyName _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  const MaterialProperty<RankTwoTensor> & _F_old;
  const std::string _ep_name;
  MaterialProperty<Real> & _ep;
  const MaterialProperty<Real> & _ep_old;
  MaterialProperty<RankTwoTensor> & _be;
  const MaterialProperty<RankTwoTensor> & _be_old;
  MaterialProperty<RankTwoTensor> & _Np;
  MaterialProperty<RankTwoTensor> &_Fp;
  const MaterialProperty<RankTwoTensor> &_Fp_old;
  MaterialProperty<RankTwoTensor> &_Fe;
  const MaterialProperty<RankTwoTensor> &_Fe_old;

  //rates

  MaterialProperty<RankTwoTensor> &_Ee;
  MaterialProperty<RankTwoTensor> &_Ee_dot;

  MaterialProperty<RankTwoTensor> &_Ep;
  MaterialProperty<RankTwoTensor> &_Ep_dot;

  MaterialBase * _flow_stress_material;
  const std::string _flow_stress_name;
  const MaterialProperty<Real> & _H;
  const MaterialProperty<Real> & _dH;
  const MaterialProperty<Real> & _d2H;

  const ADMaterialProperty<Real> &_rho;
  const Real _C0;
  const Real _C1;
  const Real _Le;
  //MaterialProperty<RankTwoTensor> &_pressure_av;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient_old;

  const MaterialProperty<RankTwoTensor> &_cauchy_stress;

  ////////////

  //fracture stuff
  const VariableValue &_c;
  const Real _l;
  MaterialProperty<Real> &_kappa;
  MaterialProperty<Real> &_L;
  const Real _gc;
  const VariableValue &_gcprop;
  const Real _visco;
  const Real _kdamage;
  MaterialProperty<Real> &_Hist;
  const MaterialProperty<Real> &_Hist_old;

  MaterialProperty<Real> &_elastic_energy;
  MaterialProperty<Real> &_delastic_energydc;
  MaterialProperty<Real> &_d2elastic_energyd2c;
  MaterialProperty<Real> &_dstress_dc;

  MaterialProperty<RankTwoTensor> &_sigma;
  MaterialProperty<Real> &_sigma_pressure;
  MaterialProperty<RankTwoTensor> &_sigma_dev;

  MaterialProperty<RankTwoTensor> &_sigma_pos;
  MaterialProperty<RankTwoTensor> &_sigma_neg;
  MaterialProperty<Real> &_W;
  MaterialProperty<Real> &_Wpos;
  MaterialProperty<Real> &_Wneg;
  //invariants for debugging
  MaterialProperty<Real> &_I1_pos;
  MaterialProperty<Real> &_I3_pos;
  MaterialProperty<Real> &_I1_neg;  
  MaterialProperty<Real> &_I3_neg;
  MaterialProperty<Real> &_elastic_energy_total;
  const Real _ep_ref;
  MaterialProperty<RankTwoTensor> &_Cp_bar;
  MaterialProperty<RankTwoTensor> &_Cp;
  const MaterialProperty<RankTwoTensor> &_Cp_bar_old;
  const MaterialProperty<RankTwoTensor> &_Cp_old;
  const MaterialProperty<RankTwoTensor> &_Ep_old;
  const MaterialProperty<RankTwoTensor> &_Ee_old;
  MaterialProperty<RankTwoTensor> &_F_computed;
  MaterialProperty<RankTwoTensor> &_S;
  MaterialProperty<Real> &_HS_elastic;
  MaterialProperty<Real> &_HS_plastic;
  MaterialProperty<RankTwoTensor> &_C_computed;
  const VariableValue &_h_min;
private:
  /// @{ Helper (dummy) variables for iteratively updating the consistant tangent during return mapping
  RankFourTensor _d_be_d_F;
  RankFourTensor _d_n_d_be;
  RankTwoTensor _d_deltaep_d_betr;
  RankTwoTensor _d_R_d_betr;
  RankTwoTensor _d_J_d_betr;
  /// @}
};