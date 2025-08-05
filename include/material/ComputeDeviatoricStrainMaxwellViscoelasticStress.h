
#pragma once

#include "ComputeStressBase.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

class ComputeDeviatoricStrainMaxwellViscoelasticStress : public ComputeStressBase
{
public:
  static InputParameters validParams();

  ComputeDeviatoricStrainMaxwellViscoelasticStress(const InputParameters & parameters);

//  void initialSetup() override;

protected:
  /// Material property defining buik modulus
  const Real _Kv;

  /// Material property defining infinity shear modulus
  const Real _uinf;

  /// Material property defining relaxation modulus
  const Real _ui;

  /// Material property defining relaxation time
  const Real _Ti;

  /// Creep strain
  MaterialProperty<RankTwoTensor> & _creep_strain;
  const MaterialProperty<RankTwoTensor> & _creep_strain_old;

  /// Other strains
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  /// Rank two identity tensor
  const RankTwoTensor _identity_two;

  /// Rank four symmetric identity tensor
  const RankFourTensor _identity_symmetric_four;

  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;
};
