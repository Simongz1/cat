#include "ADComputeMieGruneisenPressure.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeMieGruneisenPressure);

InputParameters
ADComputeMieGruneisenPressure::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes the standard form of the Mie-Gruneisen EOS pressure");

    params.addRequiredParam<Real>("K0", "name of the reference bulk modulus of the material");
    params.addRequiredParam<Real>("s", "slope of the us-up relation");
    params.addRequiredParam<Real>("gamma", "gruneisen parameter");
    params.addRequiredParam<MooseFunctorName>("density_name", "name of the density variable");
    params.addRequiredParam<MooseFunctorName>("reference_density_name", "name of the reference density variable");
    params.addRequiredParam<MooseFunctorName>("sie_name", "name of the specific internal energy");
    params.addRequiredParam<MooseFunctorName>("sie0_name", "name of the reference specific internal energy");
    params.addRequiredParam<MooseFunctorName>("pressure_name", "name to give to this property");
    params.addParam<Real>("density_limit", 1e-11 ,"minimum density value");
    return params;
}

ADComputeMieGruneisenPressure::ADComputeMieGruneisenPressure(const InputParameters &params)
    : FunctorMaterial(params),
      _K0(getParam<Real>("K0")),
      _s(getParam<Real>("s")),
      _gamma(getParam<Real>("gamma")),
      _rho(getFunctor<ADReal>("density_name")),
      _rho0(getFunctor<ADReal>("reference_density_name")),
      _sie(getFunctor<ADReal>("sie_name")),
      _sie0(getFunctor<ADReal>("sie0_name")),
      _density_limit(getParam<Real>("density_limit"))
{

    addFunctorProperty<ADReal>(
        getParam<MooseFunctorName>("pressure_name"),
        [this](const auto & r, const auto & state) -> ADReal{

            //forward declare the required quantities
            const Real K0 = _K0;
            const Real s = _s;
            const Real gamma = _gamma;
            const ADReal rho = MetaPhysicL::max(_rho(r, state), ADReal(_density_limit));
            const ADReal rho0 = _rho0(r, state);
            const ADReal sie = _sie(r, state);
            const ADReal sie0 = _sie0(r, state);
            const Real limit = _density_limit;

            //form cold curve
            ADReal x = 1. - (rho0 / rho);
            ADReal pcold = K0 * x;
            pcold *= 1. / MetaPhysicL::pow(1. - (s * x), 2.);
            pcold *= (1. - (gamma / 2.) * (rho / rho0 - 1.));

            //form thermal pressure
            ADReal phot = gamma * MetaPhysicL::max(sie - sie0, ADReal(0.));

            //add contributions

            ADReal ptot = pcold + phot;
            return ptot;
        }
    );
}