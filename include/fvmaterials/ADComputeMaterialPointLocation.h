#pragma once

#include "FunctorMaterial.h"
#include <cmath>

class ADComputeMaterialPointLocation : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeMaterialPointLocation(const InputParameters & params);
protected:
    const Moose::Functor<ADRealVectorValue> &_rhoX;
    const Moose::Functor<ADReal> &_rho;
};