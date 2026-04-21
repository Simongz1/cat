#include "ADComputeEnergyFlux.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeEnergyFlux);

InputParameters
ADComputeEnergyFlux::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes energy flux vector");
    params.addRequiredParam<MooseFunctorName>("energy", "name of the energy functor");
    params.addRequiredParam<MooseFunctorName>("density", "name of the density variable");
    params.addRequiredParam<MooseFunctorName>("mx", "x component of momentum");
    params.addRequiredParam<MooseFunctorName>("my", "y component of momentum");
    params.addRequiredParam<MooseFunctorName>("mz", "z component of momentum");
    params.addRequiredParam<MooseFunctorName>("pressure", "name of the pressure functor");
    params.addRequiredParam<MooseFunctorName>("energy_flux_name", "name of the energy flux vector");
    return params;
}

ADComputeEnergyFlux::ADComputeEnergyFlux(const InputParameters &params)
    : FunctorMaterial(params),
      _energy(getFunctor<ADReal>("energy")),
      _density(getFunctor<ADReal>("density")),

      //get velocity components
      _mx(getFunctor<ADReal>("mx")),
      _my(getFunctor<ADReal>("my")),
      _mz(getFunctor<ADReal>("mz")),
      _pressure(getFunctor<ADReal>("pressure"))
{

    //recycle this object to declare rho * u

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("energy_flux_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{
            const ADReal energy = _energy(r, state);
            const ADReal density = _density(r, state);
            const ADReal pressure = _pressure(r, state);

            const ADReal mx = _mx(r, state);
            const ADReal my = _my(r, state);
            const ADReal mz = _mz(r, state);

            //define vector valued u for easy index operation
            ADRealVectorValue m;
            m(0) = mx;
            m(1) = my;
            m(2) = mz;

            //compute the energy flux vector
            ADRealVectorValue eflux;
            for (unsigned int j = 0; j < 3; ++j){
                eflux(j) = (energy + pressure) / density;
                eflux(j) *= m(j);
            }

            return eflux;
        }
    );
}