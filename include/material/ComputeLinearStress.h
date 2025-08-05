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

class ComputeLinearStress : public ComputeStressBase
{
public:
  static InputParameters validParams();
  ComputeLinearStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress() override;
  virtual void initQpStatefulProperties() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

// temperature
  const Real _ref_temperature;
  const VariableValue & _temperature;
  const Real _alpha;
  
  const MaterialProperty<RankTwoTensor> & _total_strain;
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
  
  usingTensorIndices(i_, j_, k_, l_);
};