
#pragma once

#include "ComputeStressBase.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include "MathUtils.h"
#include "RankFourTensor.h"

class MieGruneisenStress : public ComputeStressBase
{
public:
  static InputParameters validParams();

  MieGruneisenStress(const InputParameters & parameters);

//  void initialSetup() override;

protected:
  const Real _Gamma;
  const Real _T_ref;
  const MaterialProperty<Real> &_rho;
  const MaterialProperty<Real> &_Cv;
  const VariableValue &_T;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain_old;
  const Real _C0;
  const Real _C1;

  //variable element length computation
  const Elem * const &_current_elem;
  const MaterialProperty<RankFourTensor> &_elasticity_tensor;
  const Real _viscosity_type;
  const Real _Le;
  MaterialProperty<Real> &_pressure_mg;
  MaterialProperty<Real> &_pressure_JWL;
  MaterialProperty<Real> &_pressure_total;
  MaterialProperty<Real> &_pressure_av;
  const Real _s;
  
  const Real _A;
  const Real _B;
  const Real _R1;
  const Real _R2;
  const Real _omega;
  
  const VariableValue &_Y4;
  const Real _use_JWL;
  MaterialProperty<Real> &_lambda;
  MaterialProperty<Real> &_rhorho;
  MaterialProperty<Real> &_dMG_dT;
  MaterialProperty<Real> &_dJWL_dT;

  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;
};
