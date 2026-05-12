#pragma once

#include "FunctorMaterial.h"

class ADComputeMieGruneisenPressure : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeMieGruneisenPressure(const InputParameters & params);
protected:
    const Real _K0;
    const Real _s;
    const Real _gamma;
    const Moose::Functor<ADReal> &_rho;
    const Moose::Functor<ADReal> &_rho0;
    const Moose::Functor<ADReal> &_sie;
    const Moose::Functor<ADReal> &_sie0;
    const Real _density_limit;
};