//this material computes the artificial viscosity stress
//called by other kernels/materials for their own computations

#include "ComputeAVStress.h"

registerMooseObject("catApp", ComputeAVStress);

InputParameters
ComputeAVStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Computes the Artificial Viscosity formulation as in"
                              "C.A Duarte 2021 JAP"
                              "The Computed value is called by several stress-temperature kernels");
  params.addParam<MaterialPropertyName>("density", "density", "Name of Material Property that provides the density");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient"); //C0 coefficient AV
  params.addRequiredParam<Real>("C1", "Landshoff coefficient"); //C1 coefficient AV
  params.addRequiredParam<Real>("Le", "Minumum element size");
  return params;
}

ComputeAVStress::ComputeAVStress(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)), //this is called to compute the bulk modulus
    _Gamma(getParam<Real>("Gamma")),
    _density(getMaterialProperty<Real>("density")),
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _ss_prop(declareProperty<Real>("ss_prop")), //sound speed property
    _stress_av(declareProperty<Real>("stress_av")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("Le"))
{
}

void
ComputeAVStress::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
}

void
ComputeAVStress::computeQpStress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  Real K0; //K0 is the bulk modulus, J_e is the volumetric change and J_e_dot is the volumetric change rate
  Real stress_av; //stress tensor to store the artificial viscosity stress
  K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2); //computation of the bulk modulus
  //store the computed BK value
  _bulk_modulus[_qp] = K0;

  //compute the required terms for AV 
  Real J_e = 1.0 + _mechanical_strain[_qp].trace();
  Real J_e_dot = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  Real ss_prop = std::sqrt(K0 / _density[_qp]); //compute sound speed from bulk modulus and density
  _ss_prop[_qp] = ss_prop; //store into material variable for outside use if needed (not right now)

  //define the artificial viscosity stress
  stress_av = (_C0 * _density[_qp] * J_e_dot * std::abs(J_e_dot) / std::pow(J_e, 2.0)) * std::pow(_Le, 2.0); //first term
  stress_av += (_C1 * _density[_qp] * ss_prop * (J_e_dot / J_e) * _Le); //second term
  _stress_av[_qp] = stress_av; //store value in material property
  //to convert to tensor, as the AV stress is symmetrical, and volumetric, it can be constructed as
  _stress[_qp] += stress_av * I2; //added the computed AV stress to the global stress variable
}