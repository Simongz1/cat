#include "HeatSource.h"
#include <RankTwoTensor.h>
#include "RankFourTensor.h"
#include <vector>

//forward declarte the class object to be acted upon

class MieGruneisenHS : public HeatSource
{
public:
  MieGruneisenHS(const InputParameters & parameters);
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
  const VariableValue &_Y4;
  const MaterialProperty<Real> &_dMG_dT;
  const MaterialProperty<Real> &_dJWL_dT;
  const Real _beta_av;
  const MaterialProperty<Real> &_pressure_av;
};
