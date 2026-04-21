#include "ADComputeMixturePressure.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeMixturePressure);

InputParameters
ADComputeMixturePressure::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes the gamma law mixture pressure for a two phase model");

    params.addRequiredParam<MooseFunctorName>("sie_mix", "name of the specific internal energy of the mixture");
    params.addRequiredParam<MooseFunctorName>("gamma1", "name of the gamma1 variable");
    params.addRequiredParam<MooseFunctorName>("gamma2", "name of the gamma2 variable");
    params.addRequiredParam<MooseFunctorName>("pi1", "name of the pi1 variable");
    params.addRequiredParam<MooseFunctorName>("pi2", "name of the pi2 variable");
    params.addRequiredParam<MooseFunctorName>("alpha1", "name of the alpha1 variable");
    params.addRequiredParam<MooseFunctorName>("artificial_viscosity", "name of the artificial viscosity variable");
    params.addRequiredParam<MooseFunctorName>("mixture_pressure_name", "name of the mixture pressure variable");
    return params;
}

ADComputeMixturePressure::ADComputeMixturePressure(const InputParameters &params)
    : FunctorMaterial(params),
      _sie_mix(getFunctor<ADReal>("sie_mix")),
      _gamma1(getFunctor<ADReal>("gamma1")),
      _gamma2(getFunctor<ADReal>("gamma2")),
      _pi1(getFunctor<ADReal>("pi1")),
      _pi2(getFunctor<ADReal>("pi2")),
      _alpha1(getFunctor<ADReal>("alpha1")),
      _art_vis(getFunctor<ADReal>("artificial_viscosity"))
{

    addFunctorProperty<ADReal>(
        getParam<MooseFunctorName>("mixture_pressure_name"),
        [this](const auto & r, const auto & state) -> ADReal{

            //forward declare the required quantities
            const ADReal sie_mix = _sie_mix(r, state);
            const ADReal gamma1 = _gamma1(r, state);
            const ADReal gamma2 = _gamma2(r, state);
            const ADReal pi1 = _pi1(r, state);
            const ADReal pi2 = _pi2(r, state);
            const ADReal alpha1 = _alpha1(r, state);
            const ADReal art_vis = _art_vis(r, state);

            ADReal p_num;
            p_num = (gamma1 - 1) * (gamma2 - 1) * sie_mix;
            p_num -= gamma1 * pi1 * (gamma2 - 1) * alpha1;
            p_num -= gamma2 * pi2 * (gamma1 - 1) * (1 - alpha1);

            ADReal p_den;
            p_den = (gamma2 - 1) * alpha1 + (gamma1 - 1) * (1 - alpha1);

            //define the output
            return (p_num / p_den);
        }
    );
}