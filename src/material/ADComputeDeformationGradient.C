#include "ADComputeDeformationGradient.h"

registerMooseObject("mlApp", ADComputeDeformationGradient);

InputParameters
ADComputeDeformationGradient::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADMaterial>::validParams();
  //params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("Returns the total lagrangian deformation gradient from strain and rotation increments");
  return params;
}

ADComputeDeformationGradient::ADComputeDeformationGradient(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ADMaterial>(parameters),
    //ADSingleVariableReturnMappingSolution(parameters),
    _strain_increment(getADMaterialProperty<RankTwoTensor>("strain_increment")),
    _rotation_increment(getADMaterialProperty<RankTwoTensor>("rotation_increment")),
    _F(declareADProperty<RankTwoTensor>("F")),
    _F_old(getMaterialPropertyOld<RankTwoTensor>("F")),
    _C(declareADProperty<RankTwoTensor>("C")),
    _epsilon(declareADProperty<RankTwoTensor>("epsilon"))
{}

void
ADComputeDeformationGradient::initialSetup()
{
  _F[_qp].setToIdentity();
}

void
ADComputeDeformationGradient::initQpStatefulProperties()
{
  ADMaterial::initQpStatefulProperties();
  _F[_qp].setToIdentity();
}

void
ADComputeDeformationGradient::computeQpProperties()
{ 
  ADRankTwoTensor I;
  I.setToIdentity();
  //form incremental deformation gradient
  ADRankTwoTensor dF = _rotation_increment[_qp] * (I + _strain_increment[_qp]);
  _F[_qp] = dF * _F_old[_qp];

  //declare C
  ADRankTwoTensor C = 0.5 * (_F[_qp].transpose() * _F[_qp] - I);
  _C[_qp] = C;

  //peform spectral decomposition of C to get epsilon
  ADRankTwoTensor Q;
  std::vector<ADReal> lam(3);
  
  C.symmetricEigenvaluesEigenvectors(lam, Q);

  //assemble log(C)
  ADRankTwoTensor logDiag;
  logDiag.zero();

  for(unsigned int i = 0; i < 3; ++i){
    logDiag(i,i) = MetaPhysicL::log(lam[i]);
  }

  //assemble logC
  ADRankTwoTensor logC = Q * logDiag * Q.transpose();

  //store
  _epsilon[_qp] = 0.5 * logC;
}