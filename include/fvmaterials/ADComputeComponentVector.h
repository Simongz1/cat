#pragma once

#include "FunctorMaterial.h"

class ADComputeComponentVector : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeComponentVector(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_x;
    const Moose::Functor<ADReal> &_y;
    const Moose::Functor<ADReal> &_z;
};