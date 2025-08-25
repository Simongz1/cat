//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeNeoHookeanPK2Fracture.h"

registerMooseObject("SolidMechanicsApp", ComputeNeoHookeanPK2Fracture);

InputParameters
ComputeNeoHookeanPK2Fracture::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ComputeLagrangianStressPK2>::validParams();
  params.addClassDescription("The Simo-Hughes style J2 plasticity.");
  params.addParam<std::string>("base_name", "The base name used to save the cracked stress");

  params.addRequiredCoupledVar("c", "Order parameter for damage");

  params.addParam<MaterialPropertyName>("elasticity_tensor", "elasticity_tensor", "The name of the elasticity tensor.");

  params.addParam<MaterialPropertyName>("kappa_name","kappa_op","Name of material property being created to store the interfacial parameter kappa");
  params.addParam<MaterialPropertyName>("mobility_name", "L", "Name of material property being created to store the mobility L");
  params.addRequiredParam<Real>("l","CrackLength");
  params.addRequiredParam<Real>("gc_prop","CrackEnergy");
  params.addRequiredParam<Real>("visco","ViscoParam");
  params.addParam<Real>("kdamage", 1e-9, "Stiffness of damaged matrix");

  params.addParam<MaterialPropertyName>("elastic_energy_name", "elastic_energy", "Name of material property storing the elastic energy");

  params.addParam<Real>("C0", 0, "C0 Artificial Visosity Parameter");
  params.addParam<Real>("C1", 0, "C1 Artificial Visosity Parameter");
  params.addRequiredCoupledVar("elementsize","elementsize");

  return params;
}

ComputeNeoHookeanPK2Fracture::ComputeNeoHookeanPK2Fracture( const InputParameters & parameters)
  : DerivativeMaterialInterface<ComputeLagrangianStressPK2>(parameters),
    GuaranteeConsumer(this),

    _elasticity_tensor_name(_base_name + getParam<MaterialPropertyName>("elasticity_tensor")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>(_elasticity_tensor_name)),

    _kappa(declareProperty<Real>(getParam<MaterialPropertyName>("kappa_name"))),
    _L(    declareProperty<Real>(getParam<MaterialPropertyName>("mobility_name"))),
    _l(      getParam<Real>("l")),
    _gc_prop(getParam<Real>("gc_prop")),
    _visco(  getParam<Real>("visco")),
    _kdamage(getParam<Real>("kdamage")),


    _c(coupledValue("c")),

    _F_old(            getMaterialPropertyOld<RankTwoTensor>(_base_name + "deformation_gradient")),

    _H(        declareProperty   <Real>("hist")),
    _H_old(getMaterialPropertyOld<Real>("hist")),

       _elastic_energy(    declareProperty          <Real>(getParam<MaterialPropertyName>("elastic_energy_name"))),
      _delastic_energydc(  declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("elastic_energy_name"),coupledName("c", 0))),
     _d2elastic_energyd2c(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("elastic_energy_name"), coupledName("c", 0), coupledName("c", 0))),
     _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", coupledName("c", 0))),
    //_d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    //_D(getMaterialProperty<Real>("D_name")),
    //_dDdc(getMaterialPropertyDerivative<Real>("D_name", coupledName("c", 0))),
    //_d2Dd2c(getMaterialPropertyDerivative<Real>("D_name", coupledName("c", 0), coupledName("c", 0))),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _elementsize(coupledValue("elementsize")),
    _rho(getMaterialProperty<Real>("density"))
{
}

void
ComputeNeoHookeanPK2Fracture::initialSetup()
{
  // Enforce isotropic elastic tensor
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("ComputeNeoHookeanFracture requires an isotropic elasticity tensor");
}

void
ComputeNeoHookeanPK2Fracture::initQpStatefulProperties()
{
  ComputeLagrangianStressPK2::initQpStatefulProperties();
  _H[_qp] = 0.0;
}

void
ComputeNeoHookeanPK2Fracture::computeQpPK2Stress()
{
  //usingTensorIndices(i_, j_, k_, l_, m_);
  usingTensorIndices(i_, j_, k_, l_);

  // Assign L and kappa;
  _kappa[_qp] = _gc_prop * _l;
  _L[_qp] = 1.0 / (_gc_prop * _visco);

  // Degredation Function
   Real c_factor = 1;
   if (_c[_qp] > 1){ c_factor = 0;} 
   Real D      = c_factor * std::pow((1-_c[_qp]),2) * (1.0 - _kdamage) + _kdamage;
   Real dDdc   = c_factor *     -2 * (1-_c[_qp])    * (1.0 - _kdamage);
   Real d2Dd2c = c_factor *  2                      * (1.0 - _kdamage);


  const Real E = ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  const Real G = ElasticityTensorTools::getIsotropicShearModulus (_elasticity_tensor[_qp]);
  const Real K = ElasticityTensorTools:: getIsotropicBulkModulus (_elasticity_tensor[_qp]);
  const Real lambda = G*(E-2*G)/(3*G-E);
  const Real mu     = G;
  const RankTwoTensor I       = RankTwoTensor::Identity();
  const RankTwoTensor F       = _F[_qp];
  const RankTwoTensor F_inv   = F.inverse();
  const RankTwoTensor F_T     = F.transpose();
  const RankTwoTensor F_inv_T = F.inverse().transpose();
  const RankTwoTensor C       = F_T*F; //Right Cauchy
  RankTwoTensor C_inv = (2 * _E[_qp] + RankTwoTensor::Identity()).inverse();
  //const RankTwoTensor C_inv   = F_inv*F_T; //Inverse Right Cauchy
  const Real J       = F.det();
  const Real J_dot   = ((_F[_qp].det() - _F_old[_qp].det()) / _dt);

  ///  Model follows from W = lambda / 2 * (ln J)^2 - mu * ln J + 1/2 * mu * (tr(C)- I)

  //~SaintVenat
  //RankTwoTensor pk2_stress = _elasticity_tensor[_qp]*_E[_qp];
  //Real elastic_energy      = pk2_stress.doubleContraction(_E[_qp]) / 2.0;
  //RankFourTensor Jacobian  = _elasticity_tensor[_qp];

  //~Neo Hookean
  Real elastic_energy      = (lambda/2)*std::pow(log(J),2) + (mu/2)*(C.trace()- 3) - mu*log(J);
  RankTwoTensor pk2_stress = lambda * log(J) * C_inv+ mu*(I-C_inv);
  RankFourTensor Jacobian  = lambda*(-2*log(J)* C_inv.times<i_, k_, l_, j_>(C_inv)+ C_inv.times<i_, j_, k_, l_>(C_inv)) + mu*2*C_inv.times<i_, k_, l_, j_>(C_inv);


  // Assign history variable
  if (elastic_energy > _H_old[_qp])
  {
    _H[_qp] = elastic_energy;
  }
  else
  {
    _H[_qp] = _H_old[_qp];
  }

  RankTwoTensor dSdc,dpk1_stressdc;
  _S                  [_qp] = D      * pk2_stress;
    dSdc                    = dDdc   * pk2_stress;
    dpk1_stressdc           = F * _S[_qp];
  _dstress_dc         [_qp] = dpk1_stressdc * F_T / J;
  _elastic_energy     [_qp] = D      *_H[_qp]      + _gc_prop * _c[_qp] * _c[_qp] / (2 * _l); //Actually this would be free energy = elastic + localFracture
  _delastic_energydc  [_qp] = dDdc   *_H[_qp]      + _gc_prop * _c[_qp]           / (    _l);
  _d2elastic_energyd2c[_qp] = d2Dd2c *_H[_qp]      + _gc_prop                     / (    _l);
  _C                  [_qp] = D      * Jacobian;


  //Add some artificial viscosity, so the simulation can deal
  Real ss = std::sqrt(K / _rho[_qp]);
  Real Le = _elementsize[_qp];
    RankTwoTensor I2 = RankTwoTensor::Identity();
  _S[_qp] += I2 * (_C0 * _rho[_qp] *      (J_dot * std::abs(J_dot) / std::pow(J, 2.0)) * std::pow(Le, 2.0));
  _S[_qp] += I2 * (_C1 * _rho[_qp] * ss * (J_dot                   /          J      ) *          Le      );
}

    //RankTwoTensor Cinv = (2 * _E[_qp] + RankTwoTensor::Identity()).inverse();
    //_S[_qp] = (_lambda[_qp] * log(_F[_qp].det()) - _mu[_qp]) * Cinv +
    //          _mu[_qp] * RankTwoTensor::Identity();
    //_C[_qp] =
    //    -2 * (_lambda[_qp] * log(_F[_qp].det()) - _mu[_qp]) * Cinv.times<i_, k_, l_, j_>(Cinv) +
    //    _lambda[_qp] * Cinv.times<i_, j_, k_, l_>(Cinv);













////---------------Have-A_Nice-Day---------------------------------------------------------------------------------------------------------------------------------
