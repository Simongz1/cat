#include "TimeDerivative.h"

// Forward Declarations
class Y1_dot;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class Y1_dot : public TimeDerivative
{
public:
  Y1_dot(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  std::string _base_name;
  const MaterialProperty<Real> &_r1;
  const MaterialProperty<Real> &_dr1dT;
  const unsigned int _TempId;
};