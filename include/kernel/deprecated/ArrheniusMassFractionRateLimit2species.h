/// Provides rate of reacting mass fraction decrease during chemical reaction

#pragma once

#include "DerivativeMaterialInterface.h"
#include "Kernel.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

/*
 * Arrhenius exponential law with an upper threshold of the reaction rate
 * two species are considered: specie number 1 is the kernel variable
 */

class ArrheniusMassFractionRateLimit2species : public DerivativeMaterialInterface<Kernel>
{
public:
  static InputParameters validParams();

  ArrheniusMassFractionRateLimit2species(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// optional parameter that allows multiple mechanics models to be defined
  std::string _base_name;

  /// Auxiliary variable: mass_fraction_1 (kernel variable)
  const VariableValue & _mass_fraction_1;

  /// MOOSE variable number for the mass fraction 1 variable
  const unsigned int _mass_fraction_1_var;

  /// name of mass fraction 1 variable
  VariableName _mass_fraction_1_name;

  /// Auxiliary variable: mass_fraction_2
  const VariableValue & _mass_fraction_2;

  /// MOOSE variable number for the mass fraction 2 variable
  const unsigned int _mass_fraction_2_var;

  /// name of mass fraction 2 variable
  VariableName _mass_fraction_2_name;

  /// Auxiliary variable: temperature
  const VariableValue & _temperature;

  /// MOOSE variable number for the temperature variable
  const unsigned int _temperature_var;

  /// name of temperature variable
  VariableName _temperature_name;

  /// Decrease rate of the mass fraction due to chemistry
  const MaterialProperty<Real> & _mass_fraction_rate;

  /// Derivative of the decrease rate of the mass fraction due to chemistry with respect to temperature
  const MaterialProperty<Real> & _dmass_fraction_rate_dtemperature;

  /// Derivative of the decrease rate of the mass fraction due to chemistry with respect to mass fraction 1
  const MaterialProperty<Real> & _dmass_fraction_rate_dmass_fraction_1; // diagonal Jacobian

  /// Derivative of the decrease rate of the mass fraction due to chemistry with respect to mass fraction 2
  const MaterialProperty<Real> & _dmass_fraction_rate_dmass_fraction_2; // Off diagonal Jacobian
};
