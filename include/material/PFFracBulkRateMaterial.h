/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#pragma once

#include "Material.h"
#include "Function.h"

/**
 * Phase-field fracture
 * This class obtains critical energy release rate (gc) value
 * Used by PFFRacBulkRate
 */
class PFFracBulkRateMaterial : public Material
{
public:
  static InputParameters validParams();

  PFFracBulkRateMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();
  /**
   * This function obtains the value of gc
   * Must be overidden by the user for heterogeneous gc
   */
//  virtual void getProp();

  ///Input parameter for homogeneous gc
  Real _gc;

  ///Material property where the gc values are stored
  MaterialProperty<Real> &_gc_prop;
  ///Function to specify varying gc
//  Function * _function_prop;

private:

};
