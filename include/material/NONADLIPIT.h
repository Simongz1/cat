#pragma once

#include "ComputeLagrangianStressPK1.h"
#include "GuaranteeConsumer.h"
#include "ElasticityTensorTools.h"
#include "SingleVariableReturnMappingSolution.h"
#include "Function.h"
#include "DerivativeMaterialInterface.h"

/* This class implements the Simo-Hughes style J2 plasticity */
class NONADLIPIT
  : public DerivativeMaterialInterface<ComputeLagrangianStressPK1>,
    public GuaranteeConsumer,
    public SingleVariableReturnMappingSolution
{
public:
  static InputParameters validParams();

  NONADLIPIT(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeQpPK1Stress() override;

  /// @{ The return mapping residual and derivative
  virtual Real computeReferenceResidual(const Real & effective_trial_stress,
                                        const Real & scalar) override;
  virtual Real computeResidual(const Real & effective_trial_stress, const Real & scalar) override;
  virtual Real computeDerivative(const Real & effective_trial_stress, const Real & scalar) override;
  virtual void preStep(const Real & scalar_old, const Real & residual, const Real & jacobian) override;
  virtual Real computeWvol(const Real & J);
  virtual Real computeWinv(const Real & x);
  /// @}

  const MaterialPropertyName _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  const MaterialProperty<RankTwoTensor> & _F;
  const MaterialProperty<RankTwoTensor> & _F_old;
  const std::string _ep_name;
  MaterialProperty<Real> & _ep;
  const MaterialProperty<Real> & _ep_old;
  MaterialProperty<Real> &_ep_dot;
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

  const MaterialProperty<Real> &_rho;
  const Real _C0;
  const Real _C1;
  const Real _Le;

  //fracture stuff

  const VariableValue &_c;

  const VariableValue &_gc;

  MaterialProperty<Real> &_Hist;
  const MaterialProperty<Real> &_Hist_old;

  MaterialProperty<Real> &_W0;
  MaterialProperty<Real> &_W;
  MaterialProperty<Real> &_Wpos;
  MaterialProperty<Real> &_Wneg;

  //invariants for debugging

  MaterialProperty<RankTwoTensor> &_inv_Cp;
  MaterialProperty<RankTwoTensor> &_Cp;

  const MaterialProperty<RankTwoTensor> &_inv_Cp_old;
  const MaterialProperty<RankTwoTensor> &_Cp_old;

  MaterialProperty<RankTwoTensor> &_Ce;
  const MaterialProperty<RankTwoTensor> &_Ce_old;

  const MaterialProperty<RankTwoTensor> &_Ep_old;
  const MaterialProperty<RankTwoTensor> &_Ee_old;
  MaterialProperty<RankTwoTensor> &_S;

  MaterialProperty<Real> &_HS_elastic;
  MaterialProperty<Real> &_HS_plastic;

  MaterialProperty<Real> &_D;
  MaterialProperty<RankTwoTensor> &_pk1_pos;
  MaterialProperty<RankTwoTensor> &_pk1_neg;
  const MaterialProperty<Real> &_kdamage;
  MaterialProperty<Real> &_Je;
  MaterialProperty<Real> &_Jp;
  MaterialProperty<RankTwoTensor> &_Fres;
  MaterialProperty<RankTwoTensor> &_E;
  MaterialProperty<RankTwoTensor> &_E_dot;
  MaterialProperty<RankTwoTensor> &_sigma;

private:
  /// @{ Helper (dummy) variables for iteratively updating the consistant tangent during return mapping
  RankFourTensor _d_be_d_F;
  RankFourTensor _d_n_d_be;
  RankTwoTensor _d_deltaep_d_betr;
  RankTwoTensor _d_R_d_betr;
  RankTwoTensor _d_J_d_betr;
  /// @}
};