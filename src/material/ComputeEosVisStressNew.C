/// isotropic EOS stress without phase field
/// Includes artificial viscosity and Mie Gruneisen Equation of State
/// Calculated plastic strain based on J2 plasticity

#include "ComputeEosVisStressNew.h"

registerMooseObject("TensorMechanicsApp", ComputeEosVisStressNew);

InputParameters
ComputeEosVisStressNew::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain"
                             "Considers Mie Gruneisen EOS and artificial viscosity damping");
  params.addRequiredParam<Real>("Gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addParam<MaterialPropertyName>("density", "density", "Name of Material Property that provides the density");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Name of Material Property that provides the density");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("slope_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("Le","Maximum element size");
  params.addRequiredParam<Real>("sound_speed","Speed of sound in the material");

  return params;
}

ComputeEosVisStressNew::ComputeEosVisStressNew(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),

    _Gamma(getParam<Real>("Gamma")),
    _density(getMaterialProperty<Real>("density")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _temperature(coupledValue("temperature")),
    _ref_temperature(getParam<Real>("reference_temperature")),
    _s(getParam<Real>("slope_UsUp")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _total_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "total_strain")),
    _total_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("Le")),
    _sound_speed(getParam<Real>("sound_speed")),
    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _pressure_eos(declareProperty<Real>("pressure_eos")),

    // contribution of each term to the general Pressure definition
    _peos1(declareProperty<Real>("peos1")),
    _peos2(declareProperty<Real>("peos2")),
    _deviatoric_strain(declareProperty<RankTwoTensor>("deviatoric_strain"))

{
}

void
ComputeEosVisStressNew::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
}

void
ComputeEosVisStressNew::computeQpStress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  //  Calculate "elastic" stresses due to equation of state
  // Calculate pressure from Mie Gruneisen (Menon, 2014), (Zhang, 2011)

  // compute the volumetric and deviatoric decomposition of strain
  Real K0, delta, eta, temperature, peos;
  RankTwoTensor stress_eos, stress, stress_cpl, total_strain, deviatoric_strain, volumetric_strain;
  K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  _bulk_modulus[_qp] = K0;

  delta = _mechanical_strain[_qp].trace();
  eta = - delta;

  temperature = _temperature[_qp];

  // compute eos pressure
  peos = - K0 * eta * (1.0 - (_Gamma * eta / 2.0)) / std::pow((1.0 - _s * eta), 2.0) - _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature);
  _pressure_eos[_qp] = peos;
  _peos1[_qp]= - K0 * eta * (1.0 - (_Gamma * eta / 2.0)) / std::pow((1.0 - _s * eta), 2.0);
  _peos2[_qp]= _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature);

  stress_eos = peos * I2;

  // coupling stress with volumetric strain induced stress removed
  stress_cpl = _elasticity_tensor[_qp] * (_elastic_strain_old[_qp] + _strain_increment[_qp]) - K0 * delta * I2;
  stress = stress_eos + stress_cpl;
  _stress[_qp] = stress;
  
  // decompose the strain increment into voluetric and deviatoric

  RankTwoTensor strain_increment, deviatoric_strain_increment, volumetric_strain_increment;
  strain_increment = _strain_increment[_qp];
  deviatoric_strain_increment = strain_increment.deviatoric();
  volumetric_strain_increment = strain_increment - deviatoric_strain_increment;

  // Calculate the elastic strain_increment
  _elastic_strain[_qp] = _mechanical_strain[_qp]; // the elastic strain is updated with the entire mechanical contribution
  
  //--------ComputeLinearElasticPFFractureStress::computeStrainVolDev----------
  // assume isotropic
  // assume vol-dev strain decomposition
  // Isotropic elasticity is assumed and should be enforced
  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const Real k = lambda + 2.0 * mu / LIBMESH_DIM;

  RankFourTensor I2I2 = I2.outerProduct(I2);

  // compute artificial viscosity contribution

  Real stress_av, J;
  Real total_strain_trace_rate, mechanical_strain_trace_rate;

  mechanical_strain_trace_rate = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  total_strain_trace_rate = (_total_strain[_qp].trace() - _total_strain_old[_qp].trace()) / _dt;
  
  J = 1 + _mechanical_strain[_qp].trace(); // the trace of the deviatoric strain is 1
  
  if (J < 1.0) { //this is done to ensure that only the compressive front of the wave is smoothed out
    stress_av = ( _C0 * _density[_qp] * mechanical_strain_trace_rate * std::abs(mechanical_strain_trace_rate) * std::pow(_Le,2.0) / std::pow(J,2.0) ) + ( _C1 * _density[_qp] * _sound_speed * mechanical_strain_trace_rate * _Le / J );
  }

  // stress_av = _density[_qp] * _C0 * std::pow(_Le, 2) * (std::pow(mechanical_strain_trace_rate.trace(), 2) - _C1 * _sound_speed * total_strain_rate.trace());
  _stress[_qp] += stress_av * I2; // append dissipative stress, which is still volumetric, to the total stress. The only deviatoric stress is J2 plastic stress

  // define some outputs to check behavior

  _deviatoric_strain[_qp] = _total_strain[_qp].deviatoric();
}
