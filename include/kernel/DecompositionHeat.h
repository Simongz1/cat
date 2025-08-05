#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <vector>

//forward declarte the class object to be acted upon

class DecompositionHeat : public HeatSource
{
public:
  DecompositionHeat(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  //virtual Real computeQpResudual();
  //virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  const std::vector<Real> _Qs;
  const MaterialProperty<std::vector<Real>> &_Ydot;
};
