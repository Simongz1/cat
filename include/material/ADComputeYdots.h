#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include <vector>
#include <cmath>

class ADComputeYdots;

class ADComputeYdots : public Material
{
public:
    ADComputeYdots(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties();

    const ADVariableValue & _T;

    const ADVariableValue & _Y1;
    const ADVariableValue & _Y2;
    const ADVariableValue & _Y3;
    const ADVariableValue & _Y4;

    const Real _Z1;
    const Real _Z2;
    const Real _Z3;

    const Real _E1;
    const Real _E2;
    const Real _E3;

    //have not only to be declared as AD, but the type has to be
    //passed as ADReal to be properly set up
    ADMaterialProperty<Real> &_Y1dot;
    ADMaterialProperty<Real> &_Y2dot;
    ADMaterialProperty<Real> &_Y3dot;
    ADMaterialProperty<Real> &_Y4dot;
    
    const Real _Rg;
    
    ADMaterialProperty<Real> &_r1;
    ADMaterialProperty<Real> &_r2;
    ADMaterialProperty<Real> &_r3;
};