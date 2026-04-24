#pragma once

#include "FunctorMaterial.h"

class ADComputePressureVector : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputePressureVector(const InputParameters & params);
protected:
    

    const Moose::Functor<ADReal> &_pressure;
    const unsigned int _component;
};