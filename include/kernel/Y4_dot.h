#include "TimeDerivative.h"

// Forward Declarations
class Y4_dot;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class Y4_dot : public TimeDerivative
{
public:
  Y4_dot(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  std::string _base_name;
  const MaterialProperty<Real> &_r3;
  const MaterialProperty<Real> &_dr3dT;
  const VariableValue &_Y3;
  const unsigned int _Y3Id;
  const unsigned int _TempId;
};

