//

#include "ComputeDeviatoricStrainMaxwellViscoelasticStress.h"

registerMooseObject("TensorMechanicsApp", ComputeDeviatoricStrainMaxwellViscoelasticStress);

InputParameters
ComputeDeviatoricStrainMaxwellViscoelasticStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Compute stress using viscoelasticity applied only to deviatoric part of strain");
  params.addRequiredParam<Real>("bulk_modulus", "bulk modulus of the incompressible material");
  params.addRequiredParam<Real>("infinity_shear_modulus", "Long term shear modulus of the material");
  params.addRequiredParam<Real>("relaxation_modulus", "shear modulus of the spring in the maxwell unit");
  params.addRequiredParam<Real>("relaxation_time", "characteristic time of the dashpot in the maxwell unit");
  return params;
}

ComputeDeviatoricStrainMaxwellViscoelasticStress::ComputeDeviatoricStrainMaxwellViscoelasticStress(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _Kv(getParam<Real>("bulk_modulus")),
    _uinf(getParam<Real>("infinity_shear_modulus")),
    _ui(getParam<Real>("relaxation_modulus")),
    _Ti(getParam<Real>("relaxation_time")),
    _creep_strain(declareProperty<RankTwoTensor>(isParamValid("base_name") ? _base_name + "_creep_strain" : "creep_strain")),
    _creep_strain_old(getMaterialPropertyOld<RankTwoTensor>(isParamValid("base_name") ? _base_name + "_creep_strain" : "creep_strain")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _identity_two(RankTwoTensor::initIdentity),
    _identity_symmetric_four(RankFourTensor::initIdentitySymmetricFour)
{
}

void
ComputeDeviatoricStrainMaxwellViscoelasticStress::initQpStatefulProperties()
{
  _creep_strain[_qp].zero();
}

void
ComputeDeviatoricStrainMaxwellViscoelasticStress::computeQpStress()
{
  // Distribute mechanical strain into two parts
  RankTwoTensor dev_strain = _mechanical_strain[_qp].deviatoric();
  RankTwoTensor dev_strain_old = _mechanical_strain_old[_qp].deviatoric();

  // Calculate creep strain in the dashpot
  Real hvis = 0.5 * ( 1.0 - exp(-_dt/_Ti) );
  _creep_strain[_qp] = exp(-_dt/_Ti) * _creep_strain_old[_qp] + hvis * ( dev_strain + dev_strain_old );

  // Assign value for elastic strain
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _creep_strain[_qp];

  // Calculate stress
  _stress[_qp] = _Kv * _mechanical_strain[_qp].tr() * _identity_two + 2.0 * _uinf * dev_strain + 2.0 * _ui * (dev_strain - _creep_strain[_qp]);

  // Compute dstress_dstrain
  RankFourTensor Ivol = _identity_two.outerProduct(_identity_two) / 3.0;
  RankFourTensor Idev = _identity_symmetric_four - Ivol;
  _Jacobian_mult[_qp] = _Kv * Ivol + 2.0 * _uinf * _identity_symmetric_four + 2.0 * _ui * (1.0 - hvis) * Idev;
}
