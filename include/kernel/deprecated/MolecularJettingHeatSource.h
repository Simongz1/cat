/// Calculates heat generated due to molecular jetting

#pragma once

#include "HeatSource.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// This kernel calculates the heat source term corresponding to molecular jetting

class MolecularJettingHeatSource : public HeatSource
{
public:
  static InputParameters validParams();

  MolecularJettingHeatSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  std::string _base_name;

  // damage phase field
  const VariableValue & _c;
  const bool _c_coupled;
  const unsigned int _c_var;

  const MaterialProperty<Real> & _molecular_jetting_heat;
  const MaterialProperty<Real> & _d_molecular_jetting_heat_dT;
  const MaterialProperty<Real> & _d_molecular_jetting_heat_dc;
};
