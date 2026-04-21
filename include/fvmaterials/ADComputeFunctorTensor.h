#pragma once

#include "FunctorMaterial.h"

class ADComputeFunctorTensor : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeFunctorTensor(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_density;
    const unsigned int _component;
    const Moose::Functor<ADReal> &_mx;
    const Moose::Functor<ADReal> &_my;
    const Moose::Functor<ADReal> &_mz;
    const Moose::Functor<ADReal> &_pressure;
};