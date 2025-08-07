#include "Material.h"
#include "RankTwoTensor.h"
#include <vector>
#include <cmath>

class ADComputePlasticWorkHeating;

class ADComputePlasticWorkHeating : public Material
{
public:
    ADComputePlasticWorkHeating(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties() override;

private:
    const MaterialProperty<RankTwoTensor> &_cauchy_stress;
    const Real _beta_p;
    const MaterialProperty<RankTwoTensor> &_Ep_dot;
    const ADMaterialProperty<Real> &_rho;
    const ADMaterialProperty<Real> &_cv;
    const MaterialProperty<RankTwoTensor> &_F;
    const bool _use_PK2;
    const VariableValue &_dirac_switch_react;
    const Real _thr_activation;
    ADMaterialProperty<Real> &_q_plastic;
    const bool _use_lump;
};