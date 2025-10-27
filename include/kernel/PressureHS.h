#include "HeatSource.h"
#include <RankTwoTensor.h>
#include "RankFourTensor.h"
#include <vector>

//forward declarte the class object to be acted upon

class PressureHS : public HeatSource
{
public:
  PressureHS(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  const VariableValue &_T;
  const MaterialProperty<RankFourTensor> &_elasticity_tensor;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain_old;
  const VariableValue &_Yfinal;
  const MaterialProperty<RankTwoTensor> &_dMG_dT;
  const MaterialProperty<RankTwoTensor> &_dJWL_dT;
  const Real _beta_av;
  const MaterialProperty<RankTwoTensor> &_pressure_av;
  const MaterialProperty<RankTwoTensor> &_lagrangian_strain_rate;
};
