/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State

#pragma once

#include "ComputeStressBase.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "PiecewiseBilinear.h"
#include "RankFourTensor.h"
#include "ExternalPetscSolverApp.h"
#include "petscblaslapack.h"

class ComputeLETempStrain : public ComputeStressBase
{
public:
  static InputParameters validParams();
  ComputeLETempStrain(const InputParameters & parameters);

protected:
  virtual void computeQpStress();
  virtual void initQpStatefulProperties() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  // const MaterialProperty<Real> & _pressure;
 
  // const Real _Gamma;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<Real> & _specific_heat;
  const Real _alpha;

// temperature
  const VariableValue & _temperature;
  const Real _ref_temperature;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  
  // const Real _s;

  // const MaterialProperty<RankTwoTensor> & _strain_increment;
  // const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
  // const MaterialProperty<RankTwoTensor> & _mechanical_strain;

  MaterialProperty<Real> & _new_temp;
  MaterialProperty<RankTwoTensor> & _new_stress;

  // const Real _C0;
  // const Real _C1;
  // const Real _Le;
  // const Real _sound_speed;
  // MaterialProperty<Real> & _new_temp;
  // MaterialProperty<Real> & _bulk_modulus;
  // MaterialProperty<Real> & _strain_temp;
  // MaterialProperty<Real> & _peos1;
  // MaterialProperty<Real> & _peos2;






  usingTensorIndices(i_, j_, k_, l_);
};