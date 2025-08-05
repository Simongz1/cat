/// This file includes artificial viscosity with anisortopic crack propagation

#include "ComputeLinearElasticPFFractureStressArtVisco.h"

registerMooseObject("TensorMechanicsApp", ComputeLinearElasticPFFractureStressArtVisco);


InputParameters
ComputeLinearElasticPFFractureStressArtVisco::validParams()
{
  InputParameters params = ComputeLinearElasticPFFractureStress::validParams();
  //InputParameters params = validParams<ComputeLinearElasticPFFractureStress>();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("Le","Maximum element size");
  params.addRequiredParam<Real>("sound_speed","Speed of sound in the material");
  params.addRequiredParam<Real>("density","Density of the material");
  return params;
}

ComputeLinearElasticPFFractureStressArtVisco::ComputeLinearElasticPFFractureStressArtVisco(const InputParameters & parameters)
  : ComputeLinearElasticPFFractureStress(parameters),
    _density(getParam<Real>("density")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("Le")),
    _ss(getParam<Real>("sound_speed")),
    _identity_two(RankTwoTensor::initIdentity)
{
}
void
ComputeLinearElasticPFFractureStressArtVisco::initQpStatefulProperties()
{
    //ComputePFFractureStressBase
//   ComputeLinearElasticPFFractureStress::initQpStatefulProperties();
  //_plastic_strain[_qp].zero();
  //_eqv_plastic_strain[_qp] = 0.0;

  _H[_qp] = 0.0;
  const double B=1E10;
  const double PI=3.141592653589793238463;
  // define the initial phase field history variable


 if(std::abs(_q_point[_qp](0)-12.5)<2){if(std::abs(_q_point[_qp](1)-0)<(0.3*_l[_qp])){_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* std::exp(-std::abs(_q_point[_qp](0)-12.5)/(0.3*_l[_qp])); } }

}
void
ComputeLinearElasticPFFractureStressArtVisco::computeQpStress()
{
  Real F_pos, F_neg;

  switch (_decomposition_type)
  {
    case Decomposition_type::strain_spectral:
      computeStrainSpectral(F_pos, F_neg);
      break;
    case Decomposition_type::strain_vol_dev:
      computeStrainVolDev(F_pos, F_neg);
      break;
    case Decomposition_type::stress_spectral:
      computeStressSpectral(F_pos, F_neg);
      break;
    default:
    {
      _stress[_qp] = _D[_qp] * _elasticity_tensor[_qp] * _mechanical_strain[_qp];
      F_pos = (_stress[_qp]).doubleContraction(_mechanical_strain[_qp]) / 2.0;
      F_neg = 0.0;
      if (_use_current_hist)
        _d2Fdcdstrain[_qp] = _stress[_qp] * _dDdc[_qp];

      _dstress_dc[_qp] = _stress[_qp] * _dDdc[_qp];
      _Jacobian_mult[_qp] = _D[_qp] * _elasticity_tensor[_qp];
    }
  }

  // Calculate bulk-viscosity stress term
  Real trD = _strain_increment[_qp].tr() / _dt ;
  Real q_bv = _density * _Le * ( _C0 * _Le * trD * trD - _C1 * _ss * trD );
  RankTwoTensor pressure_BV = q_bv * _identity_two;
  _stress[_qp] = _stress[_qp] - pressure_BV;

  // // Assign history variable
  if (F_pos > _H_old[_qp])
    _H[_qp] = F_pos;
  else
    _H[_qp] = _H_old[_qp];

  Real hist_variable = _H_old[_qp];
  if (_use_current_hist)
    hist_variable = _H[_qp];

  if (hist_variable < _barrier[_qp])
    hist_variable = _barrier[_qp];

  // Elastic free energy density
  _E[_qp] = hist_variable * _D[_qp] + F_neg;
  _dEdc[_qp] = hist_variable * _dDdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp];
}
