#pragma once

#include "FunctorMaterial.h"

class ADComputeEnergyFlux : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeEnergyFlux(const InputParameters & params);
protected:
    const Moose::Functor<ADReal> &_energy;
    const Moose::Functor<ADReal> &_density;

    const Moose::Functor<ADReal> &_mx;
    const Moose::Functor<ADReal> &_my;
    const Moose::Functor<ADReal> &_mz;

    const Moose::Functor<ADReal> &_pressure;
};