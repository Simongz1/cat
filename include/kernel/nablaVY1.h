#include "TimeKernel.h"

// Forward Declarations
class nablaVY1;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class nablaVY1 : public TimeKernel
{
public:
  nablaVY1(const InputParameters & parameters);
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