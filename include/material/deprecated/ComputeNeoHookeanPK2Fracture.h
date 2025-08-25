#pragma once

#include "ComputeLagrangianStressPK2.h"
#include "GuaranteeConsumer.h"
#include "ElasticityTensorTools.h"
//#include "SingleVariableReturnMappingSolution.h"
#include "Function.h"
#include "DerivativeMaterialInterface.h"
//#include "KernelBase.h"



class ComputeNeoHookeanPK2Fracture
  : public DerivativeMaterialInterface<ComputeLagrangianStressPK2>,
  //: public T,
    public GuaranteeConsumer
{
public:
  static InputParameters validParams();

  ComputeNeoHookeanPK2Fracture(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeQpPK2Stress() override;

  const MaterialPropertyName _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _L;
  Real _l;
  Real _gc_prop;
  Real _visco;
  Real _kdamage;


  const VariableValue & _c;


  const MaterialProperty<RankTwoTensor> & _F_old;

  MaterialProperty<Real> & _H;
  const MaterialProperty<Real> & _H_old;

  MaterialProperty<Real> &   _elastic_energy;
  MaterialProperty<Real> &  _delastic_energydc;
  MaterialProperty<Real> & _d2elastic_energyd2c;

  MaterialProperty<RankTwoTensor> & _dstress_dc;

  //Artificial Viscosity

  const Real _C0;
  const Real _C1;
  const   VariableValue & _elementsize;
  const MaterialProperty<Real> & _rho;
};















///////---Have-a-Nice-Day---------------------------------------------------------------------------------------------------------------------------------------