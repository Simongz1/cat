#include "ADComputeMomentum.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputeMomentum);

InputParameters
ADComputeMomentum::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes the tensor u_i u_j - sigma_ij and exposses it as a functor to be used by the FV interface");

    params.addRequiredParam<unsigned int>("component", "velocity component for which this is meant to be used");
    params.addRequiredParam<MooseFunctorName>("mx", "x component of momentum");
    params.addRequiredParam<MooseFunctorName>("my", "y component of momentum");
    params.addRequiredParam<MooseFunctorName>("mz", "z component of momentum");
    params.addRequiredParam<MooseFunctorName>("momentum_name", "name of the momentum vector");
    return params;
}

ADComputeMomentum::ADComputeMomentum(const InputParameters &params)
    : FunctorMaterial(params),

      _component(getParam<unsigned int>("component")),

      //get velocity components
      _mx(getFunctor<ADReal>("mx")),
      _my(getFunctor<ADReal>("my")),
      _mz(getFunctor<ADReal>("mz"))
{

    //recycle this object to declare rho * u

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("momentum_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{

            //forward declare the required quantities
            const ADReal mx = _mx(r, state);
            const ADReal my = _my(r, state);
            const ADReal mz = _mz(r, state);

            //define vector valued u for easy index operation
            ADRealVectorValue m;
            m(0) = mx;
            m(1) = my;
            m(2) = mz;

            //define the output
            return m;
        }
    );
}