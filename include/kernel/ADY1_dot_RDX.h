#include "TimeDerivative.h"
#include "ADKernel.h"

// Forward Declarations
class ADY1_dot_RDX;

class ADY1_dot_RDX : public ADKernel
{
public:
  ADY1_dot_RDX(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();

private:
  const ADMaterialProperty<Real> &_r1;
};