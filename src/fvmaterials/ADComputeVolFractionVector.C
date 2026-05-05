#include "ADComputeVolFractionVector.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeVolFractionVector);

InputParameters
ADComputeVolFractionVector::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("Computes the vol fraction advection field for 2 phase flow");

    params.addRequiredParam<MooseFunctorName>("density", "name of the density variable");
    params.addRequiredParam<MooseFunctorName>("mx", "x component of momentum");
    params.addRequiredParam<MooseFunctorName>("my", "y component of momentum");
    params.addRequiredParam<MooseFunctorName>("mz", "z component of momentum");
    params.addRequiredParam<MooseFunctorName>("alpha1", "vol fraction of (assumed) solid phase");
    params.addRequiredParam<MooseFunctorName>("vol_fraction_name", "name of the vol fraction advection vector field");
    params.addParam<Real>("density_limit",1e-11,"minimum density value");
    return params;
}

ADComputeVolFractionVector::ADComputeVolFractionVector(const InputParameters &params)
    : FunctorMaterial(params),
      _density(getFunctor<ADReal>("density")),

      //get velocity components
      _mx(getFunctor<ADReal>("mx")),
      _my(getFunctor<ADReal>("my")),
      _mz(getFunctor<ADReal>("mz")),
      _alpha1(getFunctor<ADReal>("alpha1")),
      _rho0(getParam<Real>("density_limit"))
{

    //recycle this object to declare rho * u

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("vol_fraction_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{

            //forward declare the required quantities
            const ADReal density = MetaPhysicL::max(_density(r, state), _rho0);
            const ADReal mx = _mx(r, state);
            const ADReal my = _my(r, state);
            const ADReal mz = _mz(r, state);
            const ADReal alpha1 = _alpha1(r, state);

            //define vector valued u for easy index operation
            ADRealVectorValue m;
            m(0) = mx;
            m(1) = my;
            m(2) = mz;

            //define the output
            return alpha1 * m / density;
        }
    );
}