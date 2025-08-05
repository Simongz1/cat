#include "ChemHeatSourceTarver.h"

registerMooseObject("beaverApp", ChemHeatSourceTarver);

InputParameters
ChemHeatSourceTarver::validParams()
{
  //InputParameters params = validParams<HeatSource>();
  InputParameters params = HeatSource::validParams();
  params.addClassDescription("HMX chemical decomposition heat source"
                             "kernel for infinitesima strain"
                             "intended to couple with MG EOS"
                             "chemical reactions based on (Tarver 1996)");
  params.addRequiredParam<Real>("Q1", "Heat of reaction for reaction 1");
  params.addRequiredParam<Real>("Q2", "Heat of reaction for reaction 2");
  params.addRequiredParam<Real>("Q3", "Heat of reaction for reaction 3");
  return params;
}

ChemHeatSourceTarver::ChemHeatSourceTarver(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _Q1(getParam<Real>("Q1")),
    _Q2(getParam<Real>("Q2")),
    _Q3(getParam<Real>("Q3")),

    _chemY1_dot(getMaterialProperty<Real>("chemY1_dot")),
    _chemY2_dot(getMaterialProperty<Real>("chemY2_dot")),
    _chemY3_dot(getMaterialProperty<Real>("chemY3_dot")),
    _chemY4_dot(getMaterialProperty<Real>("chemY4_dot")),

    //request density for heat formulation

    _density_corr(getMaterialProperty<Real>("density_corr"))
{

}

Real
ChemHeatSourceTarver::computeQpResidual()
{
  Real hs;
  hs = - (_Q1 * _chemY1_dot[_qp]) + (_Q2 * _chemY3_dot[_qp]) + (_Q3 * _chemY4_dot[_qp]);
  hs *= _density_corr[_qp]; //multiply the heat source by the actual density, corrected to account for the shock passing by

  return hs * _test[_i][_qp];
}

Real
ChemHeatSourceTarver::computeQpJacobian()
{
  Real hs;
  hs = - _Q1 * _chemY1_dot[_qp] + _Q2 * _chemY3_dot[_qp] + _Q3 * _chemY4_dot[_qp];
  
  return 0.0 * _phi[_j][_qp] * _test[_i][_qp];
}
