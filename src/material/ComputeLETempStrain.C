// isotropic plastic elastic PFF fracture
/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State
/// Calculated plastic strain based on J2 plasticity

#include "ComputeLETempStrain.h"

registerMooseObject("TensorMechanicsApp", ComputeLETempStrain);

InputParameters
ComputeLETempStrain::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Computes the temperature increase due to strain"
                             "Considers Mie Gruneisen EOS and artificial viscosity damping");
  params.addParam<MaterialPropertyName>("density", "density", "Name of Material Property that provides the density");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Name of Material Property that provides the density");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("alpha", "thermal expansion coefficient");

  return params;
}

ComputeLETempStrain::ComputeLETempStrain(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),

    _density(getMaterialProperty<Real>("density")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _alpha(getParam<Real>("alpha")),
    _temperature(coupledValue("temperature")),
    // _temperature_old(coupledOldValue("temperature")),
    _ref_temperature(getParam<Real>("reference_temperature")),

    // _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain")),
    // _elastic_strain_current(getMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    // _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    // _stress_current(getMaterialPropertyByName<RankTwoTensor>(_base_name + "flow_stress")),

    _new_temp(declareProperty<Real>("new_temperature")),
    _new_stress(declareProperty<RankTwoTensor>("stress_w_temp"))
{
}

void
ComputeLETempStrain::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();





}

void
ComputeLETempStrain::computeQpStress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  // Calculate "elastic" stresses due to equation of state
  // Calculate pressure from Mie Gruneisen (Menon, 2014), (Zhang, 2011)
  // https://en.wikipedia.org/wiki/Mie%E2%80%93Gruneisen_equation_of_state
  Real temperature, temperature_new, delta_strain, thermal_strain, mech_strain;
  RankTwoTensor stress, mech_stress, thermal_stress;
  
  delta_strain = _mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace(); //increment on mechanical strain
  temperature = _temperature[_qp]; //current temperature
  temperature_new = temperature + (1 / _alpha) * delta_strain;
  _new_temp[_qp] = temperature + temperature_new; //update temperature field with new temperature

  // thermal_stress = _elasticity_tensor[_qp] * (_alpha * _new_temp[_qp]) * I2;
  // mech_stress = _elasticity_tensor[_qp] * _mechanical_strain[_qp] * I2;

  // stress = mech_stress + thermal_stress;
  // temperature_new = (1.0 / _alpha) * (_elasticity_tensor[_qp].inverse() * )I2.doubleContraction(_elasticity_tensor[_qp] * _elastic_strain_old[_qp]); //reference bulk modulus

  // delta = _mechanical_strain[_qp].trace();
  // eta = - delta;

  // _new_stress[_qp] = stress;
  // _strain_temp[_qp] = temperature_total;
  // Calculate the elastic strain_increment
  // _elastic_strain[_qp] = _mechanical_strain[_qp];

  //--------ComputeLinearElasticPFFractureStress::computeStrainVolDev----------
  // assume isotropic
  // assume vol-dev strain decomposition
  // Isotropic elasticity is assumed and should be enforced
  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const Real k = lambda + 2.0 * mu / LIBMESH_DIM;

  RankFourTensor I2I2 = I2.outerProduct(I2);

//------------------------------------------------------------------------------
  // Calculate bulk-viscosity stress term
  // Real trD, jacob, q_bv;
  // trD = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  // jacob = 1.0 + _mechanical_strain[_qp].trace();
  // q_bv = 0.0;
  // if (jacob < 1.0) {
  //   q_bv = ( _C0 * _density[_qp] * trD * std::abs(trD) * std::pow(_Le,2.0) / std::pow(jacob,2.0) ) + ( _C1 * _density[_qp] * _sound_speed * trD * _Le / jacob );
  // }
  // _stress[_qp] += q_bv * I2; //additional term to stress

}

