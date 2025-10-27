#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <vector>

//forward declarte the class object to be acted upon

class QYdot : public HeatSource
{
public:
  QYdot(const InputParameters & parameters);
  static InputParameters validParams();
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagonalJacobian(unsigned int jvar);

private:
  const Real &_Q1;
  const Real &_Q2;
  const Real &_Q3;

  const MaterialProperty<Real> &_Y1dot;
  const MaterialProperty<Real> &_Y2dot;
  const MaterialProperty<Real> &_Y3dot;
  const MaterialProperty<Real> &_Y4dot;
  const MaterialProperty<Real> &_r1;
  const MaterialProperty<Real> &_r2;
  const MaterialProperty<Real> &_r3;

  const MaterialProperty<Real> &_dr1dT;
  const MaterialProperty<Real> &_dr2dT;
  const MaterialProperty<Real> &_dr3dT;

  const VariableValue & _Y1;
  const VariableGradient & _grad_Y1;
  const unsigned int _Y1Id;
  const VariableValue & _Y2;
  const VariableGradient & _grad_Y2;
  const unsigned int _Y2Id;
  const VariableValue & _Y3;
  const VariableGradient & _grad_Y3;
  const unsigned int _Y3Id;
  const MaterialProperty<Real> &_rho;
};
