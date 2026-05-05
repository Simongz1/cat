#include "ADComputePressureVector.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>

registerMooseObject("mlApp", ADComputePressureVector);

InputParameters
ADComputePressureVector::validParams(){
    //inherit directly the required io stream
    InputParameters params = FunctorMaterial::validParams();

    params.addClassDescription("computes the pressure times a basis vector for calculating the pressure gradient");
    params.addRequiredParam<MooseFunctorName>("pressure", "name of the pressure functor");
    params.addRequiredParam<unsigned int>("component", "component for gradient evaluation");
    params.addRequiredParam<MooseFunctorName>("pressure_vector_name", "name of the pressure time basis vector");
    //add a sign parameter to switch from positive to negative
    params.addParam<Real>("scalar_sign", 1., "sign scalar for switching from positive to negative vector definition");
    return params;
}

ADComputePressureVector::ADComputePressureVector(const InputParameters &params)
    : FunctorMaterial(params),
      _pressure(getFunctor<ADReal>("pressure")),
      _component(getParam<unsigned int>("component")),
      _sign(getParam<Real>("scalar_sign"))
{

    //recycle this object to declare rho * u

    addFunctorProperty<ADRealVectorValue>(
        getParam<MooseFunctorName>("pressure_vector_name"),
        [this](const auto & r, const auto & state) -> ADRealVectorValue{
            const ADReal pressure = _pressure(r, state);
            const unsigned int component = _component;
            const Real sign = _sign;
            
            //return based on component
            if (_component == 0){
                return {sign * pressure, 0, 0};
            }else if (_component == 1){
                return {0, sign * pressure, 0};
            }else{
                return {0, 0, sign * pressure};
            }
        }
    );
}