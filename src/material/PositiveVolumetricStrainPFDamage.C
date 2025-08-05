/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PositiveVolumetricStrainPFDamage.h"
#include "libmesh/utility.h"

registerMooseObject("TensorMechanicsApp", PositiveVolumetricStrainPFDamage);

InputParameters
PositiveVolumetricStrainPFDamage::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Phase-field fracture model with only positive volumetric energy "
                             "contribution to damage growth-isotropic elasticity and undamaged "
                             "stress under compressive strain");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<Real>("kdamage", 1e-6, "Stiffness of damaged matrix");
  return params;
}


PositiveVolumetricStrainPFDamage::PositiveVolumetricStrainPFDamage(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _c(coupledValue("c")),
    _kdamage(getParam<Real>("kdamage")),
    _G0_pos(declareProperty<Real>("G0_pos")),
    _psi(declareProperty<Real>("psi")),
    _dstress_dc(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _dG0_pos_dstrain(declareProperty<RankTwoTensor>("dG0_pos_dstrain")),
    _etens(LIBMESH_DIM),
    _epos(LIBMESH_DIM),
    _eigval(LIBMESH_DIM),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name))
{
}

void
PositiveVolumetricStrainPFDamage::computeQpStress()
{
  updateVar();
  updateJacobian();
}

void
PositiveVolumetricStrainPFDamage::updateVar()
{
  // Isotropic elasticity is assumed
  Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  Real k = lambda + (2.0 * mu / 3.0) ;
  Real c = _c[_qp];
  Real xfac = _kdamage;
  if (c < 1.0)
    xfac += Utility::pow<2>(1.0 - c);

  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

   // Tensors of outerproduct of eigen vectors
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        _etens[i](j, k) = _eigvec(j, i) * _eigvec(k, i);

  Real etr = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    etr += _eigval[i];

  Real etrpos = (std::abs(etr) + etr) / 2.0;
  Real etrneg = (std::abs(etr) - etr) / 2.0;

  RankTwoTensor stress0pos, stress0neg;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        if (i==j)
          stress0pos(i,j) = k * etrpos;
      else
          stress0pos(i,j) = 0;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      if (i==j)
          stress0neg(i,j) = k * etrneg - 2.0 * mu * (_mechanical_strain[_qp](i,j) - (etr/3.0)) ;
        else
          stress0neg(i,j) = -2.0 * mu * _mechanical_strain[_qp](i,j);

  // Damage associated with positive component of stress
  _stress[_qp] = stress0pos * xfac - stress0neg;

  //Calculation of strain energy density
  _psi[_qp] = 0.5*(_stress[_qp].doubleContraction(_mechanical_strain[_qp]));

  // Energy with positive principal strains
  _G0_pos[_qp] = ( k / 2.0 ) * Utility::pow<2>(etrpos);

  // Used in PFFracBulkRate Jacobian
  _dG0_pos_dstrain[_qp] = stress0pos;

  // Used in StressDivergencePFFracTensors Jacobian
  if (c < 1.0)
    _dstress_dc[_qp] = -stress0pos * (2.0 * (1.0 - c));
  else
    _dstress_dc[_qp].zero();
}

void
PositiveVolumetricStrainPFDamage::updateJacobian()
{
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
