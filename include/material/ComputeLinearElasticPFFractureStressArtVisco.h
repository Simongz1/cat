/// This file includes artificial viscosity with anisortopic crack propagation

#pragma once

#include "ComputeLinearElasticPFFractureStress.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"
#include "PiecewiseBilinear.h"
#include "RankFourTensor.h"

// class ComputeLinearElasticPFFractureStressArtVisco;

// InputParameters validParams<ComputeLinearElasticPFFractureStressArtVisco>();

/**
 * Phase-field fracture
 * This class computes the stress update due to artificial viscosity damping
 */
class ComputeLinearElasticPFFractureStressArtVisco : public ComputeLinearElasticPFFractureStress
{
public:
  static InputParameters validParams();
  ComputeLinearElasticPFFractureStressArtVisco(const InputParameters & parameters);


protected:
  virtual void initQpStatefulProperties();
  /// Density
  const Real _density;

  /// Other strains
  const MaterialProperty<RankTwoTensor> & _strain_increment;

  /// Damping coefficients
  const Real _C0;
  const Real _C1;

  /// Maximum element size
  const Real _Le;

  /// Material property defining speed of sound in material
  const Real _ss;

  /// Rank two identity tensor
  const RankTwoTensor _identity_two;

  virtual void computeQpStress() override;
};

