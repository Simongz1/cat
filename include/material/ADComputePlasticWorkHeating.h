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
    const ADMaterialProperty<RankTwoTensor> &_S;
    const Real _beta_p;
    const ADMaterialProperty<RankTwoTensor> &_Ep_dot;
    const ADMaterialProperty<Real> &_Je;
    const ADVariableValue &_dirac_switch_react;
    const Real _thr_activation;
    ADMaterialProperty<Real> &_q_plastic;
};