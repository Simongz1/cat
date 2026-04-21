#pragma once

#include "FunctorMaterial.h"

class ADComputeMixturePressure : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeMixturePressure(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_sie_mix;
    const Moose::Functor<ADReal> &_gamma1;
    const Moose::Functor<ADReal> &_gamma2;
    const Moose::Functor<ADReal> &_pi1;
    const Moose::Functor<ADReal> &_pi2;
    const Moose::Functor<ADReal> &_alpha1;
    const Moose::Functor<ADReal> &_art_vis;
};