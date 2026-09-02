#include "TimeDerivative.h"

// Forward Declarations
class Y2_dot_RDX;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class Y2_dot_RDX : public TimeDerivative
{
public:
  Y2_dot_RDX(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  std::string _base_name;
  const MaterialProperty<Real> &_r2;
  const MaterialProperty<Real> &_dr2dT;
  const MaterialProperty<Real> &_r1;
  const MaterialProperty<Real> &_dr1dT;
  const VariableValue &_Y1;
  const unsigned int _Y1Id;
  const unsigned int _TempId;
};

