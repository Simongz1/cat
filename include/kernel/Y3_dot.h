#include "TimeDerivative.h"

// Forward Declarations
class Y3_dot;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class Y3_dot : public TimeDerivative
{
public:
  Y3_dot(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  std::string _base_name;
  const MaterialProperty<Real> &_r3;
  const MaterialProperty<Real> &_dr3dT;
  const MaterialProperty<Real> &_r2;
  const MaterialProperty<Real> &_dr2dT;
  const VariableValue &_Y2;
  const unsigned int _Y2Id;
  const unsigned int _TempId;
};

