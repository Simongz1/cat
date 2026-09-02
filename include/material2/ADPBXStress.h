#pragma once

#include "ADComputeStressBase.h"
#include "ElasticityTensorTools.h"
#include "ADSingleVariableReturnMappingSolution.h"
#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADPBXStress
  : public DerivativeMaterialInterface<ADComputeStressBase>,
    public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADPBXStress(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeQpStress() override;

  /// @{ The return mapping residual and derivative
  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & scalar) override;
  virtual ADReal computeResidual(const ADReal & effective_trial_stress, const ADReal & scalar) override;
  virtual ADReal computeDerivative(const ADReal & effective_trial_stress, const ADReal & scalar) override;
  virtual ADReal computeJWLPressure(const Real & A, 
                                      const Real & B,
                                      const Real & R1, 
                                      const Real & R2, 
                                      const Real & omega);
  virtual ADReal computeAVPressure();

  ADMaterialProperty<RankTwoTensor> &_F;
  const MaterialProperty<RankTwoTensor> & _F_old;
  ADMaterialProperty<RankTwoTensor> &_Fhat;
  
  //binder and RDX mechanical properties
  const ADMaterialProperty<Real> &_binder_bulk;
  const ADMaterialProperty<Real> &_binder_poisson;
  const ADMaterialProperty<Real> &_binder_yield;
  ADMaterialProperty<RankTwoTensor> &_sigma_binder;

  const ADMaterialProperty<Real> &_RDX_bulk;
  const ADMaterialProperty<Real> &_RDX_poisson;
  ADMaterialProperty<RankTwoTensor> &_sigma_RDX;

  const std::string _ep_name;
  ADMaterialProperty<Real> & _ep;
  const MaterialProperty<Real> & _ep_old;
  ADMaterialProperty<Real> &_ep_dot;

  ADMaterialProperty<RankTwoTensor> & _be_bar;
  const MaterialProperty<RankTwoTensor> & _be_bar_old;
  ADMaterialProperty<RankTwoTensor> & _Np;

  ADMaterialProperty<RankTwoTensor> &_Fp;
  const MaterialProperty<RankTwoTensor> &_Fp_old;
  ADMaterialProperty<RankTwoTensor> &_Fe;
  const MaterialProperty<RankTwoTensor> &_Fe_old;

  //rates

  ADMaterialProperty<RankTwoTensor> &_Ee;
  ADMaterialProperty<RankTwoTensor> &_Ee_dot;

  ADMaterialProperty<RankTwoTensor> &_Ep;
  ADMaterialProperty<RankTwoTensor> &_Ep_dot;

  MaterialBase * _flow_stress_material;
  const std::string _flow_stress_name;
  const ADMaterialProperty<Real> & _H;
  const MaterialProperty<Real> & _dH;
  const MaterialProperty<Real> & _d2H;

  const ADMaterialProperty<Real> &_rho;
  const Real _C0;
  const Real _C1;
  const Real _Le;

  ADMaterialProperty<Real> &_p_unreacted;
  ADMaterialProperty<Real> &_p_reacted;
  ADMaterialProperty<Real> &_p_av;
  ADMaterialProperty<Real> &_p_mix;

  //mixture parameters

  const bool _use_mixture;
  const VariableValue &_fraction_csv;

  //declare properties
  ADReal _mixture_shear;
  ADReal _mixture_bulk;

  ADMaterialProperty<RankTwoTensor> &_Cp;

  const MaterialProperty<RankTwoTensor> &_Ep_old;
  const MaterialProperty<RankTwoTensor> &_Ee_old;

  //stresses 
  ADMaterialProperty<RankTwoTensor> &_S;
  ADMaterialProperty<RankTwoTensor> &_PK1;

  ADMaterialProperty<Real> &_HS_elastic;
  ADMaterialProperty<Real> &_HS_plastic;

  ADMaterialProperty<Real> &_Je;
  ADMaterialProperty<Real> &_Jp;
  ADMaterialProperty<Real> &_J;
  
  ADMaterialProperty<RankTwoTensor> &_E;
  ADMaterialProperty<RankTwoTensor> &_E_dot;

  const ADMaterialProperty<RankTwoTensor> &_strain_increment;
  const ADMaterialProperty<RankTwoTensor> &_rotation_increment;

  ADMaterialProperty<Real> &_ss;
  const ADVariableValue &_Yinitial;

  ///equation of state parameters
  const Real _A_unreacted;
  const Real _B_unreacted;
  const Real _R1_unreacted;
  const Real _R2_unreacted;
  const Real _omega_unreacted;

  const Real _A_reacted;
  const Real _B_reacted;
  const Real _R1_reacted;
  const Real _R2_reacted;
  const Real _omega_reacted;
  
  const ADVariableValue &_temperature;
  const ADMaterialProperty<Real> &_cv;
private:
    //nothing here
    //hi
};