#pragma once

#include "FunctorMaterial.h"

class ADComputeTemperature : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeTemperature(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_energy;
    const Moose::Functor<ADReal> &_density;

    const Moose::Functor<ADReal> &_mx;
    const Moose::Functor<ADReal> &_my;
    const Moose::Functor<ADReal> &_mz;

    const Moose::Functor<ADReal> &_specific_heat;
};