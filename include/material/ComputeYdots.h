#include "Material.h"
#include "RankTwoTensor.h"
#include <vector>
#include <cmath>

class ComputeYdots;

class ComputeYdots : public Material
{
public:
    ComputeYdots(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties() override;

private:
    const VariableValue & _T;
    const VariableValue & _Y1;
    const VariableValue & _Y2;
    const VariableValue & _Y3;
    const VariableValue & _Y4;

    const Real _Z1;
    const Real _Z2;
    const Real _Z3;

    const Real _E1;
    const Real _E2;
    const Real _E3;

    MaterialProperty<Real> &_Y1dot;
    MaterialProperty<Real> &_Y2dot;
    MaterialProperty<Real> &_Y3dot;
    MaterialProperty<Real> &_Y4dot;
    
    const Real _Rg;
    
    MaterialProperty<Real> &_r1;
    MaterialProperty<Real> &_dr1dT;
    MaterialProperty<Real> &_r2;
    MaterialProperty<Real> &_dr2dT;
    MaterialProperty<Real> &_r3;
    MaterialProperty<Real> &_dr3dT;
};