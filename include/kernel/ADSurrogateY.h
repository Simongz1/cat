#include "ADKernel.h"

// Forward Declarations
class ADSurrogateY;

class ADSurrogateY : public ADKernel
{
public:
  ADSurrogateY(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual ADReal computeQpResidual();

private:
  const MaterialPropertyName _surrogate_name;
  const ADMaterialProperty<Real> &_surrogate_rate;
  const bool _positive;
};