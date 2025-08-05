/// Computes the mass fraction rate for 2 species given by chemical reactions

#include "ComputeArrheniusMassFractionRateLimit2species.h"

registerMooseObject("HeatConductionApp", ComputeArrheniusMassFractionRateLimit2species);

InputParameters
ComputeArrheniusMassFractionRateLimit2species::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addParam<std::string>("base_name","Optional parameter that allows the user to define multiple "
                               "mechanics material systems on the same block, i.e. for multiple phases");
  params.addClassDescription("Arrhenius mass fraction rate for 2 species from chemistry");
  params.addRequiredCoupledVar("mass_fraction_1","Mass fraction of the first specie");
  params.addRequiredCoupledVar("mass_fraction_2","Mass fraction of the second specie");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addParam<Real>("exponential_prefactor", 0.0, "Exponential prefactor in Arrhenius");
  params.addParam<Real>("c_1", 0.0, "linearly temperature dependent Exponential prefactor in Arrhenius");
  params.addParam<Real>("exponential_coefficient", 0.0, "Exponential coefficient in Arrhenius");
  params.addParam<Real>("exponential_factor", 0.0, "Exponential factor in Arrhenius");
  params.addParam<Real>("rate_limit", 1.0e9, "Upper threshold for the reaction rate");
  params.addParam<Real>("nu_1", 0.0, "Exponential factor of the mass fraction of the first specie");
  params.addParam<Real>("nu_2", 0.0, "Exponential factor of the mass fraction of the second specie");
  //params.addParam<MaterialPropertyName>("reaction_heat", "reaction_heat","Property name of the heat of reaction");
  params.addParam<Real>("temp_ref", 6000.0, "cutt of temp for preexonential definition change");
  return params;
}

ComputeArrheniusMassFractionRateLimit2species::ComputeArrheniusMassFractionRateLimit2species(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mass_fraction_1(coupledValue("mass_fraction_1")), // kernel variable
    _mass_fraction_1_var(coupled("mass_fraction_1")),
    _mass_fraction_1_name(getVar("mass_fraction_1", 0)->name()),
    _mass_fraction_2(coupledValue("mass_fraction_2")),
    _mass_fraction_2_var(coupled("mass_fraction_2")),
    _mass_fraction_2_name(getVar("mass_fraction_2", 0)->name()),
    _temperature(coupledValue("temperature")),
    _temperature_var(coupled("temperature")),
    _temperature_name(getVar("temperature", 0)->name()),
    _exponential_prefactor(getParam<Real>("exponential_prefactor")),
    _c_1(getParam<Real>("c_1")),
    _exponential_coefficient(getParam<Real>("exponential_coefficient")),
    _exponential_factor(getParam<Real>("exponential_factor")),
    _rate_limit(getParam<Real>("rate_limit")),
    _nu_1(getParam<Real>("nu_1")),
    _nu_2(getParam<Real>("nu_2")),
    // _reaction_heat(getMaterialProperty<Real>("reaction_heat")),
    _temp_ref(getParam<Real>("temp_ref")),
    _mass_fraction_rate(declareProperty<Real>(_base_name + "mass_fraction_rate")), // residual
    _dmass_fraction_rate_dtemperature(declarePropertyDerivative<Real>(_base_name + "mass_fraction_rate", _temperature_name)), // Off diagonal Jacobian
    _dmass_fraction_rate_dmass_fraction_1(declarePropertyDerivative<Real>(_base_name + "mass_fraction_rate", _mass_fraction_1_name)), // diagonal Jacobian
    _dmass_fraction_rate_dmass_fraction_2(declarePropertyDerivative<Real>(_base_name + "mass_fraction_rate", _mass_fraction_2_name)) // Off diagonal Jacobian
{
}

void
ComputeArrheniusMassFractionRateLimit2species::computeQpProperties()
{
  Real mass_fraction_1 = std::max(_mass_fraction_1[_qp],0.0);
  Real mass_fraction_2 = std::max(_mass_fraction_2[_qp],0.0);
  mass_fraction_1 = std::min(1.0,mass_fraction_1);
  mass_fraction_2 = std::min(1.0,mass_fraction_2);

  Real exp_reaction_rate, limit;

/* exp_reaction_rate = isParamValid("_reaction_heat") ? _reaction_heat[_qp]: _exponential_prefactor;
  if ( *_reaction_heat != NULL)*/
  //exp_reaction_rate = _reaction_heat[_qp];
  // else
  exp_reaction_rate = _exponential_prefactor;
  if ( _temperature[_qp] >= _temp_ref)
  { exp_reaction_rate += _c_1 * (_temperature[_qp] - _temp_ref); }
  exp_reaction_rate *= std::exp(_exponential_coefficient - _exponential_factor / _temperature[_qp]);

  limit = 1000.0 * exp_reaction_rate;

  if ((_exponential_prefactor + _c_1 * (_temperature[_qp] - _temp_ref)) >= 0.0) {
    exp_reaction_rate = std::min(exp_reaction_rate, limit);
  }
  else {
    exp_reaction_rate = std::max(exp_reaction_rate, limit); // negative rate limit
  }

  _mass_fraction_rate[_qp] = exp_reaction_rate;

  if (_nu_1 > 0.0) { // mass fraction 1 is present in the rate equation
    _mass_fraction_rate[_qp] *= std::pow(mass_fraction_1,_nu_1);
  }

  if (_nu_2 > 0.0) { // mass fraction 2 is present in the rate equation
    _mass_fraction_rate[_qp] *= std::pow(mass_fraction_2,_nu_2);
  }

  _dmass_fraction_rate_dmass_fraction_1[_qp] = 0.0;
  _dmass_fraction_rate_dmass_fraction_2[_qp] = 0.0;
  _dmass_fraction_rate_dtemperature[_qp] = 0.0;

  if (_fe_problem.currentlyComputingJacobian())
  {
    if (exp_reaction_rate < _rate_limit) {
      _dmass_fraction_rate_dtemperature[_qp] = _mass_fraction_rate[_qp]
                                   * (_exponential_factor / (_temperature[_qp] * _temperature[_qp]));
       if ( _temperature[_qp] >= _temp_ref)

       _dmass_fraction_rate_dtemperature[_qp] += _c_1 * _mass_fraction_rate[_qp] / (_exponential_prefactor +_c_1 * (_temperature[_qp] - _temp_ref));
    }
    if (mass_fraction_1 > 1.0e-3 && _nu_1 > 0.0) {
      _dmass_fraction_rate_dmass_fraction_1[_qp] = _mass_fraction_rate[_qp] * _nu_1 / mass_fraction_1;
    }
    if (mass_fraction_2 > 1.0e-3 && _nu_2 > 0.0) {
      _dmass_fraction_rate_dmass_fraction_2[_qp] = _mass_fraction_rate[_qp] * _nu_2 / mass_fraction_2;
    }
  }
}
