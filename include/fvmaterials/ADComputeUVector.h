#pragma once

#include "FunctorMaterial.h"

class ADComputeUVector : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeUVector(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_density;
    const Moose::Functor<ADRealVectorValue> &_m_vector;
    const Real _sign;
    const Real _rho0;
};