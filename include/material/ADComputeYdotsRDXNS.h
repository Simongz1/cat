#include "Material.h"
#include "RankTwoTensor.h"
#include <vector>
#include <cmath>

class ADComputeYdotsRDXNS;

class ADComputeYdotsRDXNS : public Material
{
public:
    ADComputeYdotsRDXNS(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties() override;

private:
    const VariableValue & _T;

    const Real _Z1;
    const Real _Z2;
    const Real _E1;
    const Real _E2;
    const Real _Rg;
    ADMaterialProperty<Real> &_r1;
    ADMaterialProperty<Real> &_r2;
    ADMaterialProperty<Real> &_Q1;
    ADMaterialProperty<Real> &_Q2;
    const Real _a1;
    const Real _b1;
    const Real _a2;
    const Real _b2;
    const Real _T_trans;
    const VariableValue &_dirac_switch_react;
    const Real _switch_react;
    const Real _rate_limit;

    //chemistry variables

    const VariableValue &_Y1;
    const VariableValue &_Y2;
    const bool _use_lump;
    const ADMaterialProperty<Real> &_rho;
    const ADMaterialProperty<Real> &_cv;
    ADMaterialProperty<Real> &_q_decomposition;
    ADMaterialProperty<Real> &_Y1_dot;
    ADMaterialProperty<Real> &_Y2_dot;
    ADMaterialProperty<Real> &_Y3_dot;
};