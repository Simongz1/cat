#include "ADComputeElasticStrainEnergy.h"

registerMooseObject("mlApp", ADComputeElasticStrainEnergy);

InputParameters
ADComputeElasticStrainEnergy::validParams()
{
  InputParameters params = DerivativeMaterialInterface<ADMaterial>::validParams();
  //params += ADSingleVariableReturnMappingSolution::validParams();
  params.addClassDescription("elastic strain energy constitutive model");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addRequiredCoupledVar("c", "c");
  //variables

  //parameters
  params.addRequiredParam<Real>("lambda", "lambda");
  params.addRequiredParam<Real>("mu", "mu");

  //for switching function
  params.addRequiredParam<Real>("c1", "c1");
  params.addRequiredParam<Real>("c2", "c2");
  params.addRequiredParam<Real>("c3", "c3");
  
  return params;
}

ADComputeElasticStrainEnergy::ADComputeElasticStrainEnergy(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ADMaterial>(parameters),
    //ADSingleVariableReturnMappingSolution(parameters),
    _temperature(adCoupledValue("temperature")),
    _c(adCoupledValue("c")),
    _c_name(coupledName("c")),
    _epsilon(getADMaterialProperty<RankTwoTensor>("epsilon")),
    _epsilon_p(getADMaterialProperty<RankTwoTensor>("epsilon_p")),
    //_epsilon_c(getADMaterialProperty<RankTwoTensor>("epsilon_c")),
    _epsilon_T(getADMaterialProperty<RankTwoTensor>("epsilon_T")),
    _alpha_thermal(getADMaterialProperty<Real>("alpha_thermal")),
    _density(getADMaterialProperty<Real>("density")),
    _lambda(getParam<Real>("lambda")),
    _mu(getParam<Real>("mu")),

    _epsilon_e(declareADProperty<RankTwoTensor>("epsilon_e")),
    _We(declareADProperty<Real>("We")),
    _small_stress(declareADProperty<RankTwoTensor>("small_stress")),
    _F(getADMaterialProperty<RankTwoTensor>("F")),
    _C(getADMaterialProperty<RankTwoTensor>("C")),
    _pk1_stress(declareADProperty<RankTwoTensor>("pk1_stress")),
    _pk2_stress(declareADProperty<RankTwoTensor>("pk2_stress")),
    _stress(declareADProperty<RankTwoTensor>("stress")),
    _Hist(declareADProperty<Real>("Hist")),
    _Hist_old(getMaterialPropertyOld<Real>("Hist")),
    _Wp(getADMaterialProperty<Real>("Wp")),
    _D(getADMaterialProperty<Real>("D")),
    _D_name("D"),
    _dD(getADMaterialProperty<Real>(derivativePropertyName(_D_name, {_c_name}))),
    _dissipation(declareADProperty<Real>("dissipation")),
    _sp(getADMaterialProperty<RankTwoTensor>("sp")),
    _epsilon_p_dot(getADMaterialProperty<RankTwoTensor>("epsilon_p_dot")),
    //_epsilon_c_dot(getADMaterialProperty<RankTwoTensor>("epsilon_c_dot")),

    _c1(getParam<Real>("c1")),
    _c2(getParam<Real>("c2")),
    _c3(getParam<Real>("c3")),
    _switch(declareADProperty<Real>("switch")),

    _sM(declareADProperty<Real>("sM")),
    _nMnM(declareADProperty<RankTwoTensor>("nMnM")),
    _s_pressure(declareADProperty<Real>("s_pressure"))
{}

void
ADComputeElasticStrainEnergy::initialSetup()
{}

void
ADComputeElasticStrainEnergy::initQpStatefulProperties()
{
  ADMaterial::initQpStatefulProperties();
  _Hist[_qp] = 0.;
  _epsilon_e[_qp].zero();
}

void
ADComputeElasticStrainEnergy::computeQpProperties()
{
  ADRankTwoTensor I;
  I.setToIdentity();

  //compute all strains to subtract from total strain
  //ADRankTwoTensor epsilon_T = _alpha_thermal[_qp] * (_temperature[_qp] - 300.) * I;

  //compute elastic strain
  ADRankTwoTensor elastic;
  elastic = _epsilon[_qp] - _epsilon_p[_qp] - _epsilon_T[_qp];
  _epsilon_e[_qp] = elastic;

  //use this strain to compute strain energy 
  _We[_qp] = _lambda * MetaPhysicL::pow(_epsilon_e[_qp].trace(), 2.);
  _We[_qp] += _mu * (_epsilon_e[_qp] * _epsilon_e[_qp]).trace();
  _We[_qp] *= 1. / _density[_qp];

  //compute small stress
  _small_stress[_qp] = _lambda * _epsilon_e[_qp].trace() * I;
  _small_stress[_qp] += _mu * _epsilon_e[_qp];

  //get PK1
  //decompose C
  ADRankTwoTensor Q;
  std::vector<ADReal> lam(3);
  _C[_qp].symmetricEigenvaluesEigenvectors(lam, Q);

  //form rank 2 tensor from proncipal directions
  std::vector<ADRankTwoTensor> M(3);
  
  for (unsigned int a = 0; a < 3; ++a){
    M[a].zero();
    for (unsigned int i = 0; i < 3; ++i){
      for (unsigned j = 0; j < 3; ++j){
        M[a](i,j) = Q(i,a) * Q(j,a);
      }
    }
  }

  //now generate the rank four tensor G
  ADRankFourTensor G;
  usingTensorIndices(i_, j_, k_, l_);
  G.zero();

  for (unsigned int a = 0; a < 3; ++a){
    for (unsigned int b = 0; b < 3; ++b){
      G += M[a].times<i_, k_, j_, l_>(M[b]);
      G += M[a].times<i_, l_, j_, k_>(M[b]);
    }
  }

  //form the left
  ADRankFourTensor left;
  left.zero();
  for (unsigned int a = 0; a < 3; ++a){
    left += (1. / lam[a]) * M[a].times<i_, j_, k_, l_>(M[a]);
  }

  //theta tensor
  ADRankTwoTensor theta_tensor;
  theta_tensor.zero();
  for (unsigned int a = 0; a < 3; ++a){
    for (unsigned int b = 0; b < 3; ++b){
      theta_tensor(a, b) = (0.5 * MetaPhysicL::log(lam[a])) - 0.5 * MetaPhysicL::log(lam[b]);
      theta_tensor(a, b) *= 1. / (lam[a] - lam[b]);
    }
  }

  //right tensor
  ADRankFourTensor right;
  right.zero();
  for (unsigned int a = 0; a < 3; ++a){
    for (unsigned int b = 0; b < 3; ++b){
      if (b < a){
        right += theta_tensor(a, b) * G;
      }
    }
  }

  //total
  ADRankFourTensor PP;
  PP = left + right;

  //then we obtain PK1
  _pk1_stress[_qp] = PP * (_D[_qp] * _small_stress[_qp]);

  //use the obtained PK1 stress to compute all other required stresses
  _pk2_stress[_qp] = _F[_qp] * _pk1_stress[_qp];
  ADReal J = _F[_qp].det();
  _stress[_qp] = (1. / J) * _pk1_stress[_qp] * _F[_qp].transpose();

  //penalize all the stress measures
  _pk1_stress[_qp] *= _D[_qp];
  _pk2_stress[_qp] *= _D[_qp];
  _stress[_qp] *= _D[_qp];

  //define history variable
  ADReal driving_energy;
  driving_energy = _dD[_qp] * _We[_qp] + _Wp[_qp];

  if (driving_energy > _Hist_old[_qp]){
    _Hist[_qp] = driving_energy;
  }else{
    _Hist[_qp] = _Hist_old[_qp];
  }

  //compute dissipation
  _dissipation[_qp] = _sp[_qp].doubleContraction(_epsilon_p_dot[_qp]);

  //compute switch function for crazing and plastic flow
  //spectral decomposition of stress

  std::vector<ADReal> s_vals(3);
  ADRankTwoTensor s_dirs;
  
  _small_stress[_qp].symmetricEigenvaluesEigenvectors(s_vals, s_dirs);

  //obtain maximum principal stress
  ADReal s_max = *std::max_element(s_vals.begin(), s_vals.end());

  _sM[_qp] = s_max;

  //obtain maximum direction
  std::vector<ADReal> s_max_dir(3);

  for (unsigned int i = 0; i < 3; ++i){
    if (s_max == s_vals[0]){
      s_max_dir[i] = s_dirs(i, 0);
    }else if(s_max == s_vals[1]){
      s_max_dir[i] = s_dirs(i, 1);
    }else{
      s_max_dir[i] = s_dirs(i, 2);
    }
  }

  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      _nMnM[_qp](i,j) = s_max_dir[i] * s_max_dir[j];
    }
  }

  //compute small stress pressure
  ADReal small_p = 0.;
  for (const ADReal s : s_vals){
    small_p += s;
  }
  small_p *= ADReal(1. / 3.);
  _s_pressure[_qp] = small_p;

  //use the maximum stress for the switching function
  _switch[_qp] = s_max - (_c1 + (_c2 / small_p) + _c3 * small_p);
}