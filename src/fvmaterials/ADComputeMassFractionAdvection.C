#include "ADComputeMassFractionAdvection.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeMassFractionAdvection);

InputParameters
ADComputeMassFractionAdvection::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("Computes the mass fraction advection field for 2 phase flow");

    params.addRequiredParam<MooseFunctorName>("density", "name of the density variable");
    params.addRequiredParam<MooseFunctorName>("mx", "x component of momentum");
    params.addRequiredParam<MooseFunctorName>("my", "y component of momentum");
    params.addRequiredParam<MooseFunctorName>("mz", "z component of momentum");
    params.addRequiredParam<MooseFunctorName>("z1rho1", "mass fraction of (assumed) solid phase");
    params.addRequiredParam<MooseFunctorName>("mass_fraction_name", "name of the mass fraction advection vector field");
    return params;
}

ADComputeMassFractionAdvection::ADComputeMassFractionAdvection(const InputParameters &params)
    : FunctorMaterial(params),
      _density(getFunctor<ADReal>("density")),

      //get velocity components
      _mx(getFunctor<ADReal>("mx")),
      _my(getFunctor<ADReal>("my")),
      _mz(getFunctor<ADReal>("mz")),
      _z1(getFunctor<ADReal>("z1rho1"))
{

    //recycle this object to declare rho * u

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("mass_fraction_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{

            //forward declare the required quantities
            const ADReal density = _density(r, state);
            const ADReal mx = _mx(r, state);
            const ADReal my = _my(r, state);
            const ADReal mz = _mz(r, state);
            const ADReal z1 = _z1(r, state);

            //define vector valued u for easy index operation
            ADRealVectorValue m;
            m(0) = mx;
            m(1) = my;
            m(2) = mz;

            //define the output
            return z1 * m; //conservative form enforced
        }
    );
}