
#include "Kernel.h"

// Forward Declarations
class HMXArrheniusY4;

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)
 */
class HMXArrheniusY4 : public Kernel
{
public:
  HMXArrheniusY4(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  std::string _base_name;
  const VariableValue & _chemY3;
  const MaterialProperty<Real> & _chemY3_dot;
  const Real _Z3;
  const Real _E3;
  const Real _R_const;
  const VariableValue & _temperature;
};

