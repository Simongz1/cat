#include "ADComputeMixtureForcingTerm.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeMixtureForcingTerm);

InputParameters
ADComputeMixtureForcingTerm::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes the gamma law mixture pressure for a two phase model");

    params.addRequiredParam<MooseFunctorName>("alpha1", "name of the specific internal energy of the mixture");
    params.addRequiredParam<MooseFunctorName>("rho_mix", "mixture denisty name");
    params.addRequiredParam<MooseFunctorName>("rho1", "name of the rho1 variable");
    params.addRequiredParam<MooseFunctorName>("c1", "name of the c1 variable");
    params.addRequiredParam<MooseFunctorName>("c2", "name of the c2 variable");
    params.addRequiredParam<MooseFunctorName>("mix_momentum_vector", "name of the mixture momentum vector");
    params.addRequiredParam<MooseFunctorName>("mixture_forcing_term_name", "name of the mixture forcing term");
    return params;
}

ADComputeMixtureForcingTerm::ADComputeMixtureForcingTerm(const InputParameters &params)
    : FunctorMaterial(params),
      _alpha1(getFunctor<ADReal>("alpha1")),
      _rho_mix(getFunctor<ADReal>("rho_mix")),
      _rho1(getFunctor<ADReal>("rho1")),
      _c1(getFunctor<ADReal>("c1")),
      _c2(getFunctor<ADReal>("c2")),
      _m_mix(getFunctor<ADRealVectorValue>("mix_momentum_vector"))
{

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("mixture_forcing_term_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{

            //forward declare the required quantities
            const ADReal alpha1 = _alpha1(r, state);
            const ADReal rho_mix = _rho_mix(r, state);
            const ADReal rho1 = _rho1(r, state);
            const ADReal c1 = _c1(r, state);
            const ADReal c2 = _c2(r, state);
            const ADRealVectorValue m_mix = _m_mix(r, state);

            //compute second phase density
            ADReal rho2;
            rho2 = (rho_mix - alpha1 * rho1) / (1 - alpha1);

            ADReal num;
            num =  rho2 * c2 * c2;

            ADReal den;
            den = alpha1 * rho2 * c2 * c2;
            den += (1 - alpha1) * rho1 * c1 * c1;

            //now multiply by the momentum vector and return

            return - num / den * m_mix / rho_mix;
        }
    );
}