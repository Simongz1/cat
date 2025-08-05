#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <vector>

//forward declarte the class object to be acted upon

class IGdotRDX : public Kernel
{
public:
  IGdotRDX(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  //virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  const Real _x0;
  const Real _y0;
  const Real _a;

  const Real _x1;
  const Real _y1;
  const Real _z1;

  const Real _x2;
  const Real _y2;
  const Real _z2;

  const Real _I;
  const Real _G1;
  const Real _G2;

  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  const MaterialProperty<Real> &_pressure_total;
  const MaterialProperty<RankTwoTensor> &_mechanical_strain;

  const Real _limI;
  const Real _limG1;
  const VariableValue &_lambda;
};
