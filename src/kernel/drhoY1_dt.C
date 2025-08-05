#include "drhoY1_dt.h"

registerMooseObject("beaverApp", drhoY1_dt);

InputParameters
drhoY1_dt::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription("Y1 species time derivate term of the conservation equation"
                              "this kernel computes only the time derivative contribution"
                              "to the residual");
  params.addParam<MaterialPropertyName>("density_t", "variable storing the variable density"
                                          "different from the nominal density rho_0");
  return params;
}

drhoY1_dt::drhoY1_dt(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _density_t(getMaterialProperty<Real>("density_t")), //time and space dependent density
    _drho_dt(getMaterialProperty<Real>("drho_dt")) //partial rho_t partial t computed elsewhere
{}

Real
drhoY1_dt::computeQpResidual()
{
  //compute the weak form contribution of the term $\frac{\partial (rho Y_1)}{\partial t}$
  //the chain rule yields -> d(rho Y1) / dt = rho dY1/dt + Y1 drho/dt
  Real res1, res2; //residuals for chain rule
  res1 = _density_t[_qp] * _u_dot[_qp];
  res2 = _drho_dt[_qp] * _u[_qp];
  Real res = res1 + res2;
  return res * _test[_i][_qp];
}
