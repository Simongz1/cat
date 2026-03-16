#include "ADComputeCrazingFlowRate.h"

registerMooseObject("mlApp", ADComputeCrazingFlowRate);

InputParameters
ADComputeCrazingFlowRate::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADMaterial>::validParams();
  //params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("Crazing flow rate constitutive model");

  //variables

  //parameters
  params.addRequiredParam<Real>("craze_flow_resistance", "craze_flow_resistance");
  params.addRequiredParam<Real>("gamma_c_0", "gamma_c_0");
  params.addRequiredParam<Real>("m_craze", "m_craze");
  
  return params;
}

ADComputeCrazingFlowRate::ADComputeCrazingFlowRate(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ADMaterial>(parameters),
    //ADSingleVariableReturnMappingSolution(parameters),

    _gamma_c_0(getParam<Real>("gamma_c_0")),
    _m_craze(getParam<Real>("m_craze")),
    _craze_flow_resistance(getParam<Real>("craze_flow_resistance")),
    _sM(getADMaterialProperty<Real>("sM")),

    //declarations
    _gamma_c_dot(declareADProperty<Real>("gamma_c_dot"))
{}

void
ADComputeCrazingFlowRate::initialSetup()
{
  _gamma_c_dot[_qp] = 0.;
}

void
ADComputeCrazingFlowRate::initQpStatefulProperties()
{
  ADMaterial::initQpStatefulProperties();
  _gamma_c_dot[_qp] = 0.;
}

void
ADComputeCrazingFlowRate::computeQpProperties()
{
  _gamma_c_dot[_qp] = _gamma_c_0 * MetaPhysicL::pow((_sM[_qp] / _craze_flow_resistance), 1. / _m_craze);
}