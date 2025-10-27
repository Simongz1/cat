#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "ElasticityTensorTools.h"
#include <vector>

//forward declarte the class object to be acted upon

class VisHS : public HeatSource
{
public:
  VisHS(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  const MaterialProperty<RankTwoTensor> &_deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> &_cauchy_stress;
  const MaterialProperty<Real> &_density;
  const MaterialProperty<RankTwoTensor> &_pressure_av;
  const MaterialProperty<RankTwoTensor> &_be;
  const MaterialProperty<RankTwoTensor> &_be_old;
  const Real _beta_av;
  const Real _beta_p;
  const MaterialProperty<RankTwoTensor> &_lagrangian_strain_rate;
  const MaterialProperty<Real> &_ep_rate;
  const MaterialProperty<RankTwoTensor> &_Np;

  const MaterialProperty<RankTwoTensor> &_plastic_strain;
  const MaterialProperty<RankTwoTensor> &_plastic_strain_old;
  const MaterialProperty<RankFourTensor> &_Cijkl;
  const MaterialProperty<Real> &_alpha;
  const VariableValue &_temperature;
  const Real _shock_heat;
  const VariableValue &_c;
};
