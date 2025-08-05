
#pragma once

#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class LinearHeatSource : public HeatSource
{
public:
  
  static InputParameters validParams();
  LinearHeatSource(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  // virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  
  std::string _base_name;

  // reference temperature with zero thermal expansion
  const VariableValue & _temperature;
  // const Real _reference_temperature;
  const MaterialProperty<Real> & _specific_heat;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor; //elasticity tensor
  const Real _alpha;

  // const MaterialProperty<Real> & _vol_strain;
  
  const MaterialProperty<RankTwoTensor> & _total_strain_old;
  const MaterialProperty<RankTwoTensor> & _total_strain;
  
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
  const MaterialProperty<RankTwoTensor> & _stress;
  
  // MaterialProperty<Real> & _effective_ps,

};

