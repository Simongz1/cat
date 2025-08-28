#include "Material.h"
#include "RankTwoTensor.h"
#include <vector>
#include <cmath>

class ADComputeElasticWorkHeating;

class ADComputeElasticWorkHeating : public Material
{
public:
    ADComputeElasticWorkHeating(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties() override;

private:
    const ADVariableValue &_T;
    const Real _beta_av;
    const ADMaterialProperty<Real> &_pressure_av;
    const ADMaterialProperty<Real> &_dP_dT;
    const MaterialProperty<RankTwoTensor> &_Ee_dot;
    const ADMaterialProperty<Real> &_rho;
    const ADMaterialProperty<Real> &_cv;
    const VariableValue &_dirac_switch_react;
    const MaterialProperty<RankTwoTensor> &_Fe;
    const Real _thr_activation;
    const bool _use_PK2;
    ADMaterialProperty<Real> &_q_elastic;
    const bool _use_lump;
};