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

    const ADMaterialProperty<Real> &_rho;
    const ADMaterialProperty<Real> &_cv;
    const MaterialPropertyName _target_name;
    const ADMaterialProperty<Real> &_target;
    const MaterialPropertyName _tau_name;
    const ADMaterialProperty<Real> &_tau;
    const bool _use_PK2;
    //declaration
    ADMaterialProperty<Real> &_q_hotspot;
    const bool _use_lump;
};