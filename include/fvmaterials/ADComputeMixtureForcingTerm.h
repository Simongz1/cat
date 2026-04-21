#pragma once

#include "FunctorMaterial.h"

class ADComputeMixtureForcingTerm : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeMixtureForcingTerm(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_alpha1;
    const Moose::Functor<ADReal> &_rho_mix;
    const Moose::Functor<ADReal> &_rho1;
    const Moose::Functor<ADReal> &_c1;
    const Moose::Functor<ADReal> &_c2;
    const Moose::Functor<ADRealVectorValue> &_m_mix;
};