#pragma once

#include "FunctorMaterial.h"

class ADComputeVectorScalar : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeVectorScalar(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_scalar;
    const Moose::Functor<ADRealVectorValue> &_vector;
};