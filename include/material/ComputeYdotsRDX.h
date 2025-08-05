#include "Material.h"
#include "RankTwoTensor.h"
#include <vector>
#include <cmath>

class ComputeYdotsRDX;

class ComputeYdotsRDX : public Material
{
public:
    ComputeYdotsRDX(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties() override;

private:
    const VariableValue & _T;
    const VariableValue & _Y1;
    const VariableValue & _Y2;
    const VariableValue & _Y3;

    const Real _Z1;
    const Real _Z2;

    const Real _E1;
    const Real _E2;

    MaterialProperty<Real> &_Y1dot;
    MaterialProperty<Real> &_Y2dot;
    MaterialProperty<Real> &_Y3dot;
    
    const Real _Rg;
    const Real _MW;
    
    MaterialProperty<Real> &_r1;
    MaterialProperty<Real> &_dr1dT;
    MaterialProperty<Real> &_r2;
    MaterialProperty<Real> &_dr2dT;

    MaterialProperty<Real> &_Q1;
    MaterialProperty<Real> &_Q2;

    const Real _a1;
    const Real _b1;
    const Real _a2;
    const Real _b2;
    const Real _T_trans;
};