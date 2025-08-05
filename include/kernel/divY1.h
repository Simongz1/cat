#include "Kernel.h"

// Forward Declarations
class divY1;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class divY1 : public Kernel
{
public:
  divY1(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  //virtual Real computeQpJacobian();
  //virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  const MaterialProperty<Real> & _density_t;
  const MaterialProperty<Real> & _dvx_dx;
  const MaterialProperty<Real> & _dvx_dy;
  const MaterialProperty<Real> & _dvy_dx;
  const MaterialProperty<Real> & _dvy_dy;
  const MaterialProperty<Real> & _drho_dx;
  const MaterialProperty<Real> & _drho_dy;
  const MaterialProperty<Real> & _drho_dt_an;
  const MaterialProperty<Real> & _gradY1;
  const MaterialProperty<Real> & _gradvx;
  const MaterialProperty<Real> & _gradvy;
  const VariableValue & _vx;
  const VariableValue & _vy;
};

