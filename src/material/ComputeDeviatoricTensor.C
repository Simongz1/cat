#include "ComputeDeviatoricTensor.h"
#include "ElasticityTensorTools.h"

registerMooseObject("SolidMechanicsApp", ComputeDeviatoricTensor);

InputParameters
ComputeDeviatoricTensor::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes the strain energy density using full elasticity tensor");
  params.addRequiredParam<std::string>("property_name", "The property name to declare");
  return params;
}

ComputeDeviatoricTensor::ComputeDeviatoricTensor(const InputParameters & parameters)
  : Material(parameters),
    _property_name(getParam<std::string>("property_name")),
    _tensorFull(getMaterialPropertyByName<RankTwoTensor>(_property_name)),
    _tensorDev(declareProperty<RankTwoTensor>(_property_name+"_dev"))
{}

void
ComputeDeviatoricTensor::computeQpProperties()
{
  _tensorDev[_qp] = _tensorFull[_qp].deviatoric();
}
