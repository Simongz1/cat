#include "ADComputeFluxVectors.h"
#include "FunctorMaterial.h"
#include "SystemBase.h"
#include <cmath>
#include "MathFVUtils.h"
#include <type_traits>
#include "MooseFunctorArguments.h"
#include "FaceInfo.h"

registerMooseObject("mlApp", ADComputeFluxVectors);

InputParameters
ADComputeFluxVectors::validParams(){

    InputParameters params = FunctorMaterial::validParams();
    params.addClassDescription("computes the inverse deformation gradient times the velocity vector, and returns a specific flux direction vector");

    //request the full material point position vector
    //these can be left as the default values provided, so no need to provide then in the input file
    params.addParam<MooseFunctorName>("X00", "X00", "name of the 00 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X01", "X01", "name of the 01 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X02", "X02", "name of the 02 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X10", "X10", "name of the 10 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X11", "X11", "name of the 11 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X12", "X12", "name of the 12 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X20", "X20", "name of the 20 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X21", "X21", "name of the 21 component of the inverse deformation gradient");
    params.addParam<MooseFunctorName>("X22", "X22", "name of the 22 component of the inverse deformation gradient");

    //also provide directly momentum components and density
    params.addRequiredParam<MooseFunctorName>("momentum_vector_name", "name of the momentum vector");
    params.addRequiredParam<MooseFunctorName>("density_name", "name of the density variable");

    //request component
    params.addRequiredParam<unsigned int>("component", "the component of the flux vector to return");
    params.addRequiredParam<unsigned int>("direction", "the direction of the flux vector to return");
    params.addParam<std::string>("flux_base","b","base name for the flux vectors");
    return params;
}

ADComputeFluxVectors::ADComputeFluxVectors(const InputParameters &params)
    : FunctorMaterial(params),

    //obtain the 9 components
      _X00(getFunctor<ADReal>("X00")),
      _X01(getFunctor<ADReal>("X01")),
      _X02(getFunctor<ADReal>("X02")),
      _X10(getFunctor<ADReal>("X10")),
      _X11(getFunctor<ADReal>("X11")),
      _X12(getFunctor<ADReal>("X12")),
      _X20(getFunctor<ADReal>("X20")),
      _X21(getFunctor<ADReal>("X21")),
      _X22(getFunctor<ADReal>("X22")),

      _m(getFunctor<ADRealVectorValue>("momentum_vector_name")),
      _rho(getFunctor<ADReal>("density_name")),
      _component(getParam<unsigned int>("component")),
      _direction(getParam<unsigned int>("direction")),
      _flux_base(getParam<std::string>("flux_base"))
{   

    //for a lambda that creates stress internally
    auto computeFullb = [this](const auto & r, const auto & state) -> ADRealVectorValue{
        //obtain individual components of the inverse deformation gradient

        const ADReal X00 = _X00(r, state);
        const ADReal X01 = _X01(r, state);
        const ADReal X02 = _X02(r, state);
        const ADReal X10 = _X10(r, state);
        const ADReal X11 = _X11(r, state);
        const ADReal X12 = _X12(r, state);
        const ADReal X20 = _X20(r, state);
        const ADReal X21 = _X21(r, state);
        const ADReal X22 = _X22(r, state);

        // const ADReal rho = _rho(r, state);
        // const ADRealVectorValue rho_grad = _rho.gradient(r, state);

        //assemble the full inverse deformation gradient
        ADRankTwoTensor X;
        X(0,0) = X00; X(0,1) = X01; X(0,2) = X02;
        X(1,0) = X10; X(1,1) = X11; X(1,2) = X12;
        X(2,0) = X20; X(2,1) = X21; X(2,2) = X22;
        
        //obtain momentum and density
        const ADRealVectorValue m = _m(r, state);
        const ADReal rho = _rho(r, state);

        //form the full vector
        ADRealVectorValue b(0.0, 0.0, 0.0);
        for (unsigned int i = 0; i < 3; ++i){
            for (unsigned int j = 0; j < 3; ++j){
                b(i) += X(i,j) * m(j) / rho;
            }
        }

        //return the full vector, outside of this helper we extract specific components and directions
        return b;
    };

    //use the specified component (x,y) to extract the component of b (x) and the
    //specific direction (y) that the flux will be on
    auto formName = [this](const unsigned int component, const unsigned int direction, const std::string base) -> std::string{
        std::string name = base + "_" + std::to_string(component) + std::to_string(direction);
        return name;
    };

    //declare the flux vector
    addFunctorProperty<ADRealVectorValue>(formName(_component, _direction, _flux_base),
        [this, computeFullb](const auto & r, const auto & state) -> ADRealVectorValue{
            //obtain component and direction
            const unsigned int component = _component;
            const unsigned int direction = _direction;

            //call the b calculation
            ADRealVectorValue b = computeFullb(r, state);

            //extract component value
            ADReal value = b(component);

            //define the return vector
            ADRealVectorValue flux_vector; 
            flux_vector = value * ADRealVectorValue(1.0, 1.0, 1.0);

            for (unsigned int i = 0; i < 3; ++i){
                if (i != direction){
                    flux_vector(i) = 0.0;
                }
            }
            
            return flux_vector;
        });
}