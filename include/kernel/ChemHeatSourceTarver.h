
#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"


// Forward Declarations
class ChemHeatSourceTarver;

class ChemHeatSourceTarver : public HeatSource
{
public:
  ChemHeatSourceTarver(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  std::string _base_name;

  const Real _Q1;
  const Real _Q2;
  const Real _Q3;

  const MaterialProperty<Real> & _chemY1_dot;
  const MaterialProperty<Real> & _chemY2_dot;
  const MaterialProperty<Real> & _chemY3_dot;
  const MaterialProperty<Real> & _chemY4_dot; 
  const MaterialProperty<Real> & _density_corr;

};

