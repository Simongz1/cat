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

class ComputeEosStress : public ComputeStressBase
{
public:
  static InputParameters validParams();
  ComputeEosStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress() override;
  virtual void initQpStatefulProperties() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  const Real _Gamma;

  // const MaterialProperty<Real> & _pressure;
 
  const MaterialProperty<Real> & _density;
  const MaterialProperty<Real> & _specific_heat;

// temperature
  const VariableValue & _temperature;
  const Real _ref_temperature;


  const Real _s;

  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  MaterialProperty<Real> & _bulk_modulus;
  MaterialProperty<Real> & _pressure_eos;
  MaterialProperty<Real> & _peos1;
  MaterialProperty<Real> & _peos2;





  usingTensorIndices(i_, j_, k_, l_);
};