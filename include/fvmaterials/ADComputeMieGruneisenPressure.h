#pragma once

#include "FunctorMaterial.h"

class ADComputeMieGruneisenPressure : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeMieGruneisenPressure(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_sie_mix;
    const Moose::Functor<ADReal> &_gamma1;
    const Moose::Functor<ADReal> &_gamma2;
    const Moose::Functor<ADReal> &_pi1;
    const Moose::Functor<ADReal> &_pi2;
    const Moose::Functor<ADReal> &_alpha1;
};