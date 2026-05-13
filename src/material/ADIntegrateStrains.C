#include "ADIntegrateStrains.h"

registerMooseObject("mlApp", ADIntegrateStrains);

InputParameters
ADIntegrateStrains::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADMaterial>::validParams();
  //params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("Euler integration of strains");
  
  //variables
  params.addRequiredCoupledVar("temperature", "temperature");

  //parameters
  return params;
}

ADIntegrateStrains::ADIntegrateStrains(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ADMaterial>(parameters),
    //ADSingleVariableReturnMappingSolution(parameters),
    _epsilon_p(declareADProperty<RankTwoTensor>("epsilon_p")),
    _epsilon_c(declareADProperty<RankTwoTensor>("epsilon_c")),
    _epsilon_T(declareADProperty<RankTwoTensor>("epsilon_T")),

    _epsilon_p_old(getMaterialPropertyOld<RankTwoTensor>("epsilon_p")),
    _epsilon_c_old(getMaterialPropertyOld<RankTwoTensor>("epsilon_c")),

    //get old here for integration
    _epsilon_p_dot(getADMaterialProperty<RankTwoTensor>("epsilon_p_dot")),
    _epsilon_c_dot(getADMaterialProperty<RankTwoTensor>("epsilon_c_dot")),
    _temperature(adCoupledValue("temperature")),
    _alpha_thermal(getADMaterialProperty<Real>("alpha_thermal"))
{}

void
ADIntegrateStrains::initialSetup()
{}

void
ADIntegrateStrains::initQpStatefulProperties()
{
  ADMaterial::initQpStatefulProperties();
  _epsilon_p[_qp].zero();
  _epsilon_c[_qp].zero();
  _epsilon_T[_qp].zero();
}

void
ADIntegrateStrains::computeQpProperties()
{
  ADRankTwoTensor I;
  I.setToIdentity();

  _epsilon_p[_qp] = _epsilon_p_old[_qp] + _epsilon_p_dot[_qp] * _dt;
  _epsilon_c[_qp] = _epsilon_c_old[_qp] + _epsilon_c_dot[_qp] * _dt;
  _epsilon_T[_qp] = _alpha_thermal[_qp] * (_temperature[_qp] - 300.) * I;
}