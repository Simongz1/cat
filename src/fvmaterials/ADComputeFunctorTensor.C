#include "ADComputeFunctorTensor.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>


registerMooseObject("mlApp", ADComputeFunctorTensor);

InputParameters
ADComputeFunctorTensor::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes the tensor u_i u_j - sigma_ij and exposses it as a functor to be used by the FV interface");

    params.addRequiredParam<MooseFunctorName>("density", "name of the density variable");
    params.addRequiredParam<unsigned int>("component", "velocity component for which this is meant to be used");
    params.addRequiredParam<MooseFunctorName>("mx", "mx momentum component");
    params.addRequiredParam<MooseFunctorName>("my", "my momentum component");
    params.addRequiredParam<MooseFunctorName>("mz", "mz momentum component");
    params.addRequiredParam<MooseFunctorName>("flux_name", "name of flux vector to generate");
    return params;
}

ADComputeFunctorTensor::ADComputeFunctorTensor(const InputParameters &params)
    : FunctorMaterial(params),
      _density(getFunctor<ADReal>("density")),
      _component(getParam<unsigned int>("component")),

      //get momentum components
      _mx(getFunctor<ADReal>("mx")),
      _my(getFunctor<ADReal>("my")),
      _mz(getFunctor<ADReal>("mz"))
{

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("flux_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{

            //store values at actual locations
            const ADReal density = _density(r, state);
            const ADReal mx = _mx(r, state);
            const ADReal my = _my(r, state);
            const ADReal mz = _mz(r, state);

            //define vector valued m for easy index operation
            ADRealVectorValue m;
            m(0) = mx;
            m(1) = my;
            m(2) = mz;

            //define outer product container
            ADRealVectorValue outer;
            outer.zero();

            //move through indexes with i fixed as the provided component
            for (unsigned int j = 0; j < 3; ++j){
                outer(j) = (1. / density) * m(_component) * m(j);
            }
            return outer;
        }
    );
}