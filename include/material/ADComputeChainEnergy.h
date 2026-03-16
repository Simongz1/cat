#pragma once

#include "ADComputeStressBase.h"
#include "ElasticityTensorTools.h"
#include "ADSingleVariableReturnMappingSolution.h"
#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class ADComputeChainEnergy
  : public DerivativeMaterialInterface<ADComputeStressBase>,
    public ADSingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  ADComputeChainEnergy(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeQpStress() override;

  /// @{ The return mapping residual and derivative
  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & scalar) override;
  virtual ADReal computeResidual(const ADReal & effective_trial_stress, const ADReal & scalar) override;
  virtual ADReal computeDerivative(const ADReal & effective_trial_stress, const ADReal & scalar) override;
  /// @}

  const MaterialPropertyName _elasticity_tensor_name;
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;

  ADMaterialProperty<RankTwoTensor> & _F;
  const MaterialProperty<RankTwoTensor> & _F_old;
  ADMaterialProperty<RankTwoTensor> &_Fhat;

  const std::string _ep_name;
  ADMaterialProperty<Real> & _ep;
  const MaterialProperty<Real> & _ep_old;
  ADMaterialProperty<Real> &_ep_dot;
  ADMaterialProperty<RankTwoTensor> & _be;
  const MaterialProperty<RankTwoTensor> & _be_old;
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

  //fracture stuff

  const ADVariableValue &_c;

  const ADVariableValue &_gc;

  ADMaterialProperty<Real> &_Hist;
  const MaterialProperty<Real> &_Hist_old;

  ADMaterialProperty<Real> &_W0;
  ADMaterialProperty<Real> &_W;
  ADMaterialProperty<Real> &_Wpos;
  ADMaterialProperty<Real> &_Wneg;

  //invariants for debugging

  ADMaterialProperty<RankTwoTensor> &_inv_Cp;
  ADMaterialProperty<RankTwoTensor> &_Cp;

  const MaterialProperty<RankTwoTensor> &_inv_Cp_old;
  const MaterialProperty<RankTwoTensor> &_Cp_old;

  ADMaterialProperty<RankTwoTensor> &_Ce;
  const MaterialProperty<RankTwoTensor> &_Ce_old;

  const MaterialProperty<RankTwoTensor> &_Ep_old;
  const MaterialProperty<RankTwoTensor> &_Ee_old;
  ADMaterialProperty<RankTwoTensor> &_S;

  ADMaterialProperty<Real> &_HS_elastic;
  ADMaterialProperty<Real> &_HS_plastic;

  const ADMaterialProperty<Real> &_D;

  const ADMaterialProperty<Real> &_kdamage;
  ADMaterialProperty<Real> &_Je;
  ADMaterialProperty<Real> &_Jp;
  ADMaterialProperty<Real> &_J;

  ADMaterialProperty<RankTwoTensor> &_Fres;
  ADMaterialProperty<RankTwoTensor> &_E;
  ADMaterialProperty<RankTwoTensor> &_E_dot;
  ADMaterialProperty<RankTwoTensor> &_sigma;

  const unsigned int _ndisp;

  const ADMaterialProperty<RankTwoTensor> &_strain_increment;
  const ADMaterialProperty<RankTwoTensor> &_rotation_increment;

  ADMaterialProperty<RankTwoTensor> &_F_from_increment;
  const bool _use_custom;
  const ADMaterialProperty<Real> &_nu;
  const ADVariableValue &_h_min;
  ADMaterialProperty<Real> &_wp;
  const MaterialProperty<Real> &_wp_old;
  ADMaterialProperty<Real> &_wp_dot;
  const ADMaterialProperty<Real> &_alpha;
  const ADVariableValue &_temperature;
  const ADMaterialProperty<Real> &_beta_heat;

  std::vector<const ADVariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_old;
};