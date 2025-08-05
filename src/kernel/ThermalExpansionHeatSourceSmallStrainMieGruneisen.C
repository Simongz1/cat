/// Calculates heat generated due to thermal expansion

#include "ThermalExpansionHeatSourceSmallStrainMieGruneisen.h"

registerMooseObject("TensorMechanicsApp",ThermalExpansionHeatSourceSmallStrainMieGruneisen);

InputParameters
ThermalExpansionHeatSourceSmallStrainMieGruneisen::validParams()
{
  InputParameters params = HeatSource::validParams();
  params.addClassDescription("Thermal expansion heat source kernel generic kernel for small strain"
                             "Mie Gruneisen equation of state (Menon, 2014) (Zhang, 2011)");
  params.addRequiredParam<Real>("Gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>("density", "density", "Property name of the density material property");
  params.addCoupledVar("c","Phase field damage variable: Used to indicate calculation of Off Diagonal Jacobian term");
  return params;
}

ThermalExpansionHeatSourceSmallStrainMieGruneisen::ThermalExpansionHeatSourceSmallStrainMieGruneisen(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _Gamma(getParam<Real>("Gamma")), // Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density")),
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _c(coupledValue("c")),
    _c_coupled(isCoupled("c")),
    _c_var(_c_coupled ? coupled("c") : 0),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _shock_heat(getMaterialProperty<Real>("shock_heat"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
ThermalExpansionHeatSourceSmallStrainMieGruneisen::computeQpResidual()
{
  return - _shock_heat[_qp] * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceSmallStrainMieGruneisen::computeQpJacobian()
{
  RankTwoTensor ee_rate = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  // approximate (small strain) Jacobian, first order in the strain rate and temperature
  // so that only the EOS term contribution has to be considered
  return _Gamma * _density[_qp] * _specific_heat[_qp] * ee_rate.trace() * _phi[_j][_qp] * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceSmallStrainMieGruneisen::computeQpOffDiagJacobian(unsigned int jvar)
{
  RankTwoTensor ee_rate = (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]) / _dt;
  Real val = 0.0;
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k]) {
      // approximate (small strain) off-diagonal Jacobian, first order in the strain rate and temperature
      val = _Gamma * _density[_qp] * _specific_heat[_qp] * _u[_qp] * _grad_phi[_j][_qp](k) * _test[_i][_qp] / _dt;
    }
  }

  if (jvar == _c_var) {
    // approximate (small strain) off-diagonal Jacobian, first order in the strain rate and temperature
    // assumption that damage is always present
    val = - 2.0 * (1.0 - _c[_qp]) * _Gamma * _density[_qp] * _specific_heat[_qp]
        * _u[_qp] * ee_rate.trace() * _phi[_j][_qp] * _test[_i][_qp];
  }

  return val;
}
