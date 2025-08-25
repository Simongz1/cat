#include "TimeKernel.h"

// Forward Declarations
class nablaVY4;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class nablaVY4 : public TimeKernel
{
public:
  nablaVY4(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();

private:
  const MaterialProperty<Real> & _rho;
  const VariableValue & _vx;
  const VariableValue & _vy;
  //gradients
  const VariableGradient & _grad_vx;
  const VariableGradient & _grad_vy;
};