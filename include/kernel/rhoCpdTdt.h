#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <vector>
#include "DerivativeMaterialInterface.h"
#include "TimeKernel.h"

//forward declarte the class object to be acted upon

class rhoCpdTdt : public DerivativeMaterialInterface<TimeKernel>
{
public:
  rhoCpdTdt(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  //virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  const MaterialProperty<Real> &_Cp;
  const MaterialProperty<Real> &_density_t;
  const std::string &_u_name;
  const MaterialProperty<Real> &_dCpdT;
  const MaterialProperty<Real> &_drho_dt;
};
