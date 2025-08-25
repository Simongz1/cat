#include "TimeDerivative.h"

// Forward Declarations
class drhoY1_dt;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class drhoY1_dt : public TimeDerivative
{
public:
  drhoY1_dt(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  //virtual Real computeQpJacobian();
  //virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  std::string _base_name;
  const MaterialProperty<Real> & _density_t;
  const MaterialProperty<Real> & _drho_dt;
};

