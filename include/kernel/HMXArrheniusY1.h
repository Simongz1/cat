
#include "Kernel.h"

// Forward Declarations
class HMXArrheniusY1;

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)
 */
class HMXArrheniusY1 : public Kernel
{
public:
  HMXArrheniusY1(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  std::string _base_name;
  const VariableValue & _chemY1;
  const MaterialProperty<Real> & _chemY1_dot;
  const Real _Z1;
  const Real _E1;
  const Real _R_const;
  const VariableValue & _temperature;
};

