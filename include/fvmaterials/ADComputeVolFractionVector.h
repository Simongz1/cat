#pragma once

#include "FunctorMaterial.h"

class ADComputeVolFractionVector : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeVolFractionVector(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_density;
    const Moose::Functor<ADReal> &_mx;
    const Moose::Functor<ADReal> &_my;
    const Moose::Functor<ADReal> &_mz;
    const Moose::Functor<ADReal> &_alpha1;
    const Real _rho0;
};