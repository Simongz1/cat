#pragma once

#include "FunctorMaterial.h"
#include <cmath>

class ADComputeFluxVectors : public FunctorMaterial
{
public:
    static InputParameters validParams();
    ADComputeFluxVectors(const InputParameters & params);
protected:
    // const Moose::Functor<ADReal> &_rhoX_x;
    // const Moose::Functor<ADReal> &_rhoX_y;
    // const Moose::Functor<ADReal> &_rhoX_z;
    const Moose::Functor<ADReal> &_X00;
    const Moose::Functor<ADReal> &_X01;
    const Moose::Functor<ADReal> &_X02;
    const Moose::Functor<ADReal> &_X10;
    const Moose::Functor<ADReal> &_X11;
    const Moose::Functor<ADReal> &_X12;
    const Moose::Functor<ADReal> &_X20;
    const Moose::Functor<ADReal> &_X21;
    const Moose::Functor<ADReal> &_X22;

    const Moose::Functor<ADRealVectorValue> &_m;
    const Moose::Functor<ADReal> &_rho;
    const unsigned int _component;
    const unsigned int _direction;
    const std::string _flux_base;
};