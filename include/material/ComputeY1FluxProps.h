#include "Material.h"
#include "RankTwoTensor.h"

class ComputeY1FluxProps;

class ComputeY1FluxProps : public Material
{
public:
    ComputeY1FluxProps(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties() override;

private:
    const MaterialProperty<Real> & _density_0;
    const MaterialProperty<RankTwoTensor> & _mechanical_strain;
    const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
    MaterialProperty<RankTwoTensor> & _epsilon_mech_dot;
    const VariableValue & _vx;
    const VariableValue & _vy;
    const VariableValue & _Y1;
    MaterialProperty<Real> & _density_t;
    const MaterialProperty<Real> & _density_t_old;
    MaterialProperty<Real> & _dvx_dx;
    MaterialProperty<Real> & _dvx_dy;
    MaterialProperty<Real> & _dvy_dx;
    MaterialProperty<Real> & _dvy_dy;
    MaterialProperty<Real> & _drho_dx;
    MaterialProperty<Real> & _drho_dy;
    MaterialProperty<Real> & _drho_dt;
    MaterialProperty<Real> & _drho_dt_an;
    const VariableGradient & _gradY1;
    MaterialProperty<Real> & _dY1_dx;
    MaterialProperty<Real> & _dY1_dy;
    const VariableGradient & _gradvx;
    const VariableGradient & _gradvy;
};