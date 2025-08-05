#include "ADKernel.h"

// Forward Declarations
class ADY2_dot_RDX;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class ADY2_dot_RDX : public ADKernel
{
public:
  ADY2_dot_RDX(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();

private:
  const ADMaterialProperty<Real> &_r2;
  const ADMaterialProperty<Real> &_r1;
  const VariableValue &_Y1;
};

