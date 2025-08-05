/// isotropic plastic elastic PFF fracture
/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State
/// Calculated plastic strain based on J2 plasticity

#include "ComputeLinearStress.h"


registerMooseObject("TensorMechanicsApp", ComputeLinearStress);

InputParameters
ComputeLinearStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Computes linear elastic stress for incremental small strain formulation, includes thermal eigenstrain contribution"
                             "Considers isotropic elasticity tensor");
  params.addRequiredParam<Real>("alpha", "thermal_expansion_coefficient");
  params.addRequiredParam<Real>("reference_temperature", "reference_temperature");
  // params.addRequiredParam<Real>("Le","Maximum element size");
  // params.addRequiredParam<Real>("sound_speed","Speed of sound in the material");
  params.addRequiredCoupledVar("temperature", "temperature");
  return params;
}

ComputeLinearStress::ComputeLinearStress(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),

    _ref_temperature(getParam<Real>("reference_temperature")),
    _temperature(coupledValue("temperature")),
    _alpha(getParam<Real>("alpha")),
    
    _total_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain")),    
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain"))

{
}

void
ComputeLinearStress::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
}

void
ComputeLinearStress::computeQpStress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  //  Calculate "elastic" stresses due to equation of state
  // Calculate pressure from Mie Gruneisen (Menon, 2014), (Zhang, 2011)
  RankTwoTensor stress, sb;
  RankFourTensor C;

  C = _elasticity_tensor[_qp];
  sb = _total_strain[_qp] - ((_alpha) * (_temperature[_qp] - _ref_temperature)) * I2;
  stress = C * sb;

  // update stress on qp

  _stress[_qp] = stress;
 
  // Calculate the elastic strain_increment
  _elastic_strain[_qp] = _mechanical_strain[_qp];
}
