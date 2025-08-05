/// Calculates heat generated due to thermal expansion

#include "MolecularJettingHeatSource.h"

registerMooseObject("TensorMechanicsApp",MolecularJettingHeatSource);

InputParameters
MolecularJettingHeatSource::validParams()
{
  InputParameters params = HeatSource::validParams();
  params.addCoupledVar("c","Phase field damage variable");
  params.addClassDescription("Molecular jetting hear source kernel");
  return params;
}

MolecularJettingHeatSource::MolecularJettingHeatSource(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _c(coupledValue("c")),
    _c_coupled(isCoupled("c")),
    _c_var(_c_coupled ? coupled("c") : 0),
    _molecular_jetting_heat(getMaterialProperty<Real>("molecular_jetting_heat")),
    _d_molecular_jetting_heat_dT(getMaterialProperty<Real>("d_molecular_jetting_heat_dT")),
    _d_molecular_jetting_heat_dc(getMaterialProperty<Real>("d_molecular_jetting_heat_dc"))
{
}

Real
MolecularJettingHeatSource::computeQpResidual()
{
  return - _molecular_jetting_heat[_qp] * _test[_i][_qp];
}

Real
MolecularJettingHeatSource::computeQpJacobian()
{
  return - _d_molecular_jetting_heat_dT[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}

Real
MolecularJettingHeatSource::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val = 0.0;
  if (jvar == _c_var) {
    val = _d_molecular_jetting_heat_dc[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  }
  return val;
}
