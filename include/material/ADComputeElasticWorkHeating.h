#include "ADMaterial.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"
#include <vector>
#include <cmath>

class ADComputeElasticWorkHeating : public DerivativeMaterialInterface<Material>
{
public:
    ADComputeElasticWorkHeating(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties() override;

private:
    const ADVariableValue &_T;
    const ADVariableGradient &_Tgrad;
    const Real _beta_av;
    const ADMaterialProperty<RankTwoTensor> &_S;
    const ADMaterialProperty<Real> &_dPdT;
    const ADMaterialProperty<RankTwoTensor> &_Ee_dot;
    const ADMaterialProperty<Real> &_rho;
    const ADMaterialProperty<Real> &_cv;
    const ADVariableValue &_dirac_switch_react;
    const ADMaterialProperty<RankTwoTensor> &_Fe;
    const Real _thr_activation;

    ADMaterialProperty<Real> &_q_elastic;
    ADMaterialProperty<Real> &_norm_gradT;

};