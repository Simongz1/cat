#pragma once

#include "FunctorMaterial.h"
#include <cmath>

class ADComputeStressTensor : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeStressTensor(const InputParameters & params);
protected:
    const Moose::Functor<ADRealVectorValue> &_grad_x;
    const Moose::Functor<ADRealVectorValue> &_grad_y;
    const Moose::Functor<ADRealVectorValue> &_grad_z;
    const Moose::Functor<ADReal> &_lambda;
    const Moose::Functor<ADReal> &_mu;
    const Moose::Functor<ADRealVectorValue> &_m;
    const Moose::Functor<ADReal> &_rho;
};