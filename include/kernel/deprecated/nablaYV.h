#include "TimeKernel.h"

// Forward Declarations
class nablaYV;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class nablaYV : public TimeKernel
{
public:
  nablaYV(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const MaterialProperty<Real> & _rho;
  const VariableValue & _vx;
  const VariableValue & _vy;
  //gradients
  const VariableGradient & _grad_vx;
  const VariableGradient & _grad_vy;
};