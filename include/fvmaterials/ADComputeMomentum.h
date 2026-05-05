#pragma once

#include "FunctorMaterial.h"

class ADComputeMomentum : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeMomentum(const InputParameters & params);
protected:
    const unsigned int _component;
    const Moose::Functor<ADReal> &_mx;
    const Moose::Functor<ADReal> &_my;
    const Moose::Functor<ADReal> &_mz;
};