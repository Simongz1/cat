/// isotropic plastic elastic PFF fracture
/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State
/// Calculated plastic strain based on J2 plasticity

#include "ComputeEosPFStress.h"

registerMooseObject("TensorMechanicsApp", ComputeEosPFStress);

InputParameters
ComputeEosPFStress::validParams()
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
  params.addRequiredParam<Real>("k_dmg", "residual stiffness for c = 1");
  params.addRequiredParam<Real>("Le","Maximum element size");
  params.addRequiredParam<Real>("ss","Speed of sound in the material");
  params.addRequiredCoupledVar("c", "c");

  params.addParam<MaterialPropertyName>("E_name", "elastic_energy", "name of elastic energy material prop");
  params.addParam<MaterialPropertyName>("F_name", "local_fracture_energy", "name of local fracture energy material prop");
  params.addParam<MaterialPropertyName>("D_name", "degradation", "name of degradation material prop");
  params.addParam<MaterialPropertyName>("I_name", "indicator", "name of indicator material prop");
  params.addParam<bool>("use_current_history_value", false, "use the current value for the history variable");

  return params;
}

ComputeEosPFStress::ComputeEosPFStress(const InputParameters & parameters)
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
    _total_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain")),
    _total_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "total_strain")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("Le")),
    _ss(getParam<Real>("ss")),
    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _pressure_eos(declareProperty<Real>("pressure_eos")),

    // contribution of each term to the general Pressure definition
    _peos1(declareProperty<Real>("peos1")),
    _peos2(declareProperty<Real>("peos2")),
    _c(coupledValue("c")),
    _k_dmg(getParam<Real>("k_dmg")),

    // damage stress parameters F, D, I, E

    _E(declareProperty<Real>(getParam<MaterialPropertyName>("E_name"))),
    _dEdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name())),
    _d2Ed2c(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _I(getDefaultMaterialProperty<Real>("I_name")),
    _dIdc(getMaterialPropertyDerivative<Real>("I_name", getVar("c", 0)->name())),
    _d2Id2c(getMaterialPropertyDerivative<Real>("I_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _D(getMaterialProperty<Real>("D_name")),
    _dDdc(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name())),
    _d2Dd2c(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _f_pressure(getDefaultMaterialProperty<Real>("fracture_pressure")),
    _H(declareProperty<Real>("hist")),
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _use_current_hist(getParam<bool>("use_current_history_value")),
    _elastic_energy(declareProperty<Real>("elastic_energy")),
    _deviatoric_strain(declareProperty<RankTwoTensor>("deviatoric_strain"))

{
}

void
ComputeEosPFStress::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
}

void
ComputeEosPFStress::computeQpStress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  //  Calculate "elastic" stresses due to equation of state
  // Calculate pressure from Mie Gruneisen (Menon, 2014), (Zhang, 2011)
  Real K0, delta, eta, temperature, peos;
  RankTwoTensor stress_eos, stress, stress_cpl;
  K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  _bulk_modulus[_qp] = K0;
  delta = _total_strain[_qp].trace(); // use total strain as it considers volumetric expansion effects (which are promoted by temperature)
  eta = - delta;
  temperature = _temperature[_qp];
  peos = - K0 * eta * (1.0 - (_Gamma * eta / 2.0)) / std::pow((1.0 - _s * eta), 2.0) - _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature);
  _pressure_eos[_qp] = peos;
  _peos1[_qp] = - K0 * eta * (1.0 - (_Gamma * eta / 2.0)) / std::pow((1.0 - _s * eta), 2.0);
  _peos2[_qp] = _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature);

  stress_eos = peos * I2;

  // coupling stress with volumetric strain induced stress removed

  stress_cpl = _elasticity_tensor[_qp] * (_elastic_strain_old[_qp] + _strain_increment[_qp]) - K0 * delta * I2;
  stress = stress_eos + stress_cpl; // until this point stress is only eos minus the volumetric contribution

  //--------ComputeLinearElasticPFFractureStress::computeStrainVolDev----------
  // assume isotropic
  // assume vol-dev strain decomposition
  // Isotropic elasticity is assumed and should be enforced
  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const Real k = lambda + 2.0 * mu / LIBMESH_DIM;

  RankFourTensor I2I2 = I2.outerProduct(I2);

  // compute artificial viscosity contribution

  Real stress_av, J, epsilon_m_trace_dot;
  RankTwoTensor epsilon_T_dot;
  
  J = 1 + _mechanical_strain[_qp].trace();
  epsilon_m_trace_dot = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt; // use total strain as temperature can generate oscilations on the stress profile

  stress_av = stress_av = ( _C0 * _density[_qp] * epsilon_m_trace_dot * std::abs(epsilon_m_trace_dot) * std::pow(_Le,2.0) / std::pow(J,2.0) ) + ( _C1 * _density[_qp] * _ss * epsilon_m_trace_dot * _Le / J );
  
  _stress[_qp] = stress;
  
  if (J < 1.0) {
      _stress[_qp] += std::abs(stress_av) * I2;
  }
  
  _elastic_strain[_qp] = _total_strain[_qp] - _total_strain[_qp].deviatoric();

  // NOTE
  // the av stress is not yet added to the total stress as first it
  // will be pensalized by the degradation function

  //___________________

  // from this part the PF definition starts
  // first, we define the spectral decomposition for energy, stress and strain

  RankTwoTensor stress0, stress0pos, stress0neg, elastic_strain_dev, elastic_strain_vol;
  Real elastic_strain_trace, elastic_strain_trace_pos, elastic_strain_trace_neg;
  Real F_pos, F_neg;

  elastic_strain_trace = _elastic_strain[_qp].trace();

  elastic_strain_trace_pos = (std::abs(elastic_strain_trace) + elastic_strain_trace) / 2;
  elastic_strain_trace_neg = (std::abs(elastic_strain_trace) - elastic_strain_trace) / 2;

  elastic_strain_dev = _elastic_strain[_qp].deviatoric(); // computes deviatoric part of elastic strain
  elastic_strain_vol = _elastic_strain[_qp] - elastic_strain_dev; // computes volumetric elastic strain component

  stress0pos = K0 * elastic_strain_trace_pos * I2 + mu * elastic_strain_dev; //compressive component of stress (which is eos + cpl + av)
  stress0neg = stress - stress0pos; // here stress has only Peos + Pcpl, still has not been modified with av, should also have plasticity from J2 kernels

  F_pos = 0.5 * K0 * std::pow(elastic_strain_trace_pos, 2) + mu * elastic_strain_dev.doubleContraction(elastic_strain_dev);
  F_neg = 0.5 * K0 * std::pow(elastic_strain_trace_neg, 2);

  RankTwoTensor stress_orig;
  
  stress_orig = _stress[_qp]; // take the original stress before penalization

  _stress[_qp] = stress0pos * _D[_qp] - _f_pressure[_qp] * I2 * _I[_qp] + stress0neg; //penalizes the tensile component of stress with updated value of c damage variable
  _stress[_qp] += stress_av * I2; // add to the updated stress the artificial viscosity contribution

  // deffinition of history variable value

  Real hist_var = _H_old[_qp];

  if (F_pos > _H_old[_qp])
    _H[_qp] = F_pos;
  else
    _H[_qp] = _H_old[_qp];
  
  if (_use_current_hist)
    hist_var = _H[_qp];
  
  // compute energy derivatives and strain energy as they are requested by ComputePFFractureMechanicsOffDiag modules

  _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp]; // second derivative of fracture energy wrt c and epsilon
  _E[_qp] = hist_var * _D[_qp] + F_neg - _f_pressure[_qp] * _elastic_strain[_qp].trace() * _I[_qp];
  _dEdc[_qp] = hist_var * _dDdc[_qp] - _f_pressure[_qp] * _elastic_strain[_qp].trace() * _dIdc[_qp];
  _d2Ed2c[_qp] = hist_var * _d2Dd2c[_qp] - _f_pressure[_qp] * _elastic_strain[_qp].trace() * _d2Id2c[_qp];
  _elastic_energy[_qp] = _E[_qp]; // compute elastic energy requested by fracture driving energy
  _deviatoric_strain[_qp] = (_elastic_strain[_qp]).deviatoric();
}
