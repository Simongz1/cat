#include "Material.h"
#include "RankTwoTensor.h"
#include <vector>
#include <cmath>

class ADComputeHotspotHeatRate;

class ADComputeHotspotHeatRate : public Material
{
public:
    ADComputeHotspotHeatRate(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties() override;

private:
    const VariableValue &_T;
    const ADMaterialProperty<Real> &_rho;
    const ADMaterialProperty<Real> &_cv;
    const VariableValue &_target;
    const VariableValue &_tau_hotspot;
    const Real _T_ref;
    //declaration
    ADMaterialProperty<Real> &_q_hotspot;
    const bool _use_lump;
    const bool _direct_T;
    const bool _use_current;
    const Real _k_proportional;
    const VariableValue &_dirac_switch_react;
    const bool _use_sin;
    const Real _pi;
};