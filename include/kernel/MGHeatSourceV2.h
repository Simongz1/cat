#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

//forward declarte the class object to be acted upon
class MGHeatSourceV2;

class MGHeatSourceV2 : public HeatSource
{
public:
  MGHeatSourceV2(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  //virtual Real computeQpResudual();
  //virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  
  std::string _base_name;

  // Gruneisen G (or Gamma) parameter (Menon, 2014)
  const Real _gamma;
  // reference temperature with zero thermal expansion

  const MaterialProperty<Real> & _specific_heat;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor; //elasticity tensor
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
  const MaterialProperty<RankTwoTensor> & _elastic_strain;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const Real _beta_av;
  const MaterialProperty<Real> & _stress_av;
};
