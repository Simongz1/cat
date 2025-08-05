#include "ADKernel.h"

// Forward Declarations
class ADY3_dot_RDX;

//this kernel computes the term $\frac{\partial (rho Y_1)}{\partial t}$
//term for the species conservation equation

class ADY3_dot_RDX : public ADKernel
{
public:
  ADY3_dot_RDX(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();

private:
  const ADMaterialProperty<Real> &_r2;
  const VariableValue &_Y2;
};

