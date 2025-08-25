/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State

#include "ComputeStressBase.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "MathUtils.h"
#include "petscblaslapack.h"


class ComputeAVStress : public ComputeStressBase
{
public:
  ComputeAVStress(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual void computeQpStress() override;
  virtual void initQpStatefulProperties() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  const Real _Gamma;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  //declared Mat.Prop
  MaterialProperty<Real> & _bulk_modulus;
  MaterialProperty<Real> & _ss_prop;
  MaterialProperty<Real> & _stress_av;
  const Real _C0;
  const Real _C1;
  const Real _Le;
  usingTensorIndices(i_, j_, k_, l_);
};