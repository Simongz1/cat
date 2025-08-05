
#include "Kernel.h"

// Forward Declarations
class HMXArrheniusY2;

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)
 */
class HMXArrheniusY2 : public Kernel
{
public:
  HMXArrheniusY2(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  std::string _base_name;
  const VariableValue & _chemY2;
  const MaterialProperty<Real> & _chemY2_dot;
  const Real _Z2;
  const Real _E2;
  const Real _R_const;
  const VariableValue & _temperature;
  const Real _Z1;
  const Real _E1;
  const VariableValue & _chemY1;
};

