#include "ADComputeTemperature.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeTemperature);

InputParameters
ADComputeTemperature::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes energy flux vector");
    params.addRequiredParam<MooseFunctorName>("energy", "name of the energy functor");
    params.addRequiredParam<MooseFunctorName>("density", "name of the density variable");
    params.addRequiredParam<MooseFunctorName>("mx", "x component of momentum");
    params.addRequiredParam<MooseFunctorName>("my", "y component of momentum");
    params.addRequiredParam<MooseFunctorName>("mz", "z component of momentum");
    params.addRequiredParam<MooseFunctorName>("specific_heat", "name of the energy flux vector");
    params.addRequiredParam<MooseFunctorName>("specific_internal_energy", "name of the energy flux vector");
    params.addRequiredParam<MooseFunctorName>("temperature", "name of the energy flux vector");
    
    return params;
}

ADComputeTemperature::ADComputeTemperature(const InputParameters &params)
    : FunctorMaterial(params),
      _energy(getFunctor<ADReal>("energy")), //this provides rho * E
      _density(getFunctor<ADReal>("density")),

      //get velocity components
      _mx(getFunctor<ADReal>("mx")),
      _my(getFunctor<ADReal>("my")),
      _mz(getFunctor<ADReal>("mz")),
      _specific_heat(getFunctor<ADReal>("specific_heat"))
{

    //recycle this object to declare rho * u

    addFunctorProperty<ADReal>(
        getParam<MooseFunctorName>("specific_internal_energy"),
        [this](const auto & r, const auto & state) -> ADReal{
            const ADReal energy = _energy(r, state);
            const ADReal density = _density(r, state);

            const ADReal mx = _mx(r, state);
            const ADReal my = _my(r, state);
            const ADReal mz = _mz(r, state);

            //define vector valued u for easy index operation
            ADRealVectorValue m;
            m(0) = mx;
            m(1) = my;
            m(2) = mz;

            //compute specific internal energy
            ADReal sp_int_e;
            sp_int_e = energy;
            for (unsigned int i = 0; i < 3; ++i){
                sp_int_e -= (density / 2.) * (m(i) / density) * (m(i) / density);
            }
            sp_int_e *= 1. / density;

            return sp_int_e;
        }
    );

    addFunctorProperty<ADReal>(
        getParam<MooseFunctorName>("temperature"),
        [this](const auto & r, const auto & state) -> ADReal{
            const ADReal energy = _energy(r, state);
            const ADReal density = _density(r, state);

            const ADReal mx = _mx(r, state);
            const ADReal my = _my(r, state);
            const ADReal mz = _mz(r, state);
            const ADReal cv = _specific_heat(r, state);

            //define vector valued u for easy index operation
            ADRealVectorValue m;
            m(0) = mx;
            m(1) = my;
            m(2) = mz;

            //compute specific internal energy
            ADReal sp_int_e;
            sp_int_e = energy;
            for (unsigned int i = 0; i < 3; ++i){
                sp_int_e -= (density / 2.) * (m(i) / density) * (m(i) / density);
            }
            sp_int_e *= 1. / density;

            //use a linear relation between u and T with cv
            ADReal T;
            T = sp_int_e / cv;
            return T;
        }
    );
}