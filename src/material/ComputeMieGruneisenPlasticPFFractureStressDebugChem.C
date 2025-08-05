/// Computes EoS stress with artificial viscosity and J2 plasticity with Johnson - Cook model for yield. Radial Return algorithm for returning to the yield surface
/// SGZ 2024

#include "ComputeMieGruneisenPlasticPFFractureStressDebugChem.h"

registerMooseObject("TensorMechanicsApp", ComputeMieGruneisenPlasticPFFractureStressDebugChem);

InputParameters
ComputeMieGruneisenPlasticPFFractureStressDebugChem::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain"
                             "Considers Mie Gruneisen EOS and artificial viscosity damping");
  params.addRequiredParam<Real>("Gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addParam<MaterialPropertyName>("density", "density", "Name of Material Property that provides the density");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Name of Material Property that provides the density");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("slope_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("Le","Maximum element size");
  params.addParam<Real>("rtol", 1e-8, "Plastic strain NR tolerance");
  params.addParam<Real>("ftol", 1e-4, "Consistency condition NR tolerance");
  params.addParam<Real>("eptol", 1e-7, "Equivalent plastic strain NR tolerance");
  params.addRequiredParam<Real>("A_yield", "A coefficient for yield stress");
  params.addRequiredParam<Real>("B_yield", "B coefficient for yield stress");
  params.addRequiredParam<Real>("C_yield", "C coefficient for yield stress");
  params.addRequiredParam<Real>("a_corr_yield", "small a correction for yield stress");
  params.addRequiredParam<Real>("m_yield", "m coefficient for temperature for yield stress");
  params.addRequiredParam<Real>("n_yield", "n coefficient for temperature for yield stress");
  params.addRequiredParam<Real>("strain_ref", "reference plastic strain rate for yield stress");
  params.addRequiredParam<Real>("T_melt_ref", "reference melting temperature for yield stress");
  params.addRequiredParam<Real>("maxiter", "number of maximum iterations for the Radial Return algorithm convergence");
  params.addRequiredParam<Real>("toly", "tolerance for convergence on Radial Return algorithm");
  params.addRequiredParam<Real>("k_alpha", "factor for modifying the undamaged alpha parameter");
  params.addRequiredParam<Real>("thr_gas", "damage threshold for gas compression stress source activation");
  params.addRequiredParam<Real>("thr_linear", "thr_linear");
  params.addRequiredParam<Real>("yield3_thr", "threshold for thermal softening");
  params.addRequiredParam<Real>("gamma_pd", "phono drag limit value");

  //chemical reactions params

  params.addRequiredCoupledVar("chemY1", "mass fraction Y1");
  params.addRequiredCoupledVar("chemY2", "mass fraction Y2");
  params.addRequiredCoupledVar("chemY3", "mass fraction Y3");
  params.addRequiredCoupledVar("chemY4", "mass fraction Y4");

  params.addRequiredParam<Real>("Z1", "Z1 chemical parameter");
  params.addRequiredParam<Real>("Z2", "Z2 chemical parameter");
  params.addRequiredParam<Real>("Z3", "Z3 chemical parameter");
  params.addRequiredParam<Real>("Z4", "Z4 chemical parameter");

  params.addRequiredParam<Real>("E1", "E1 chemical parameter");
  params.addRequiredParam<Real>("E2", "E2 chemical parameter");
  params.addRequiredParam<Real>("E3", "E3 chemical parameter");

  params.addRequiredParam<Real>("R_const", "gas constant");

  //parameters for JWL EOS for the gas

  params.addRequiredParam<Real>("A_g", "A_g parameter for JWL eos");
  params.addRequiredParam<Real>("B_g", "B_g parameter for JWL eos");
  params.addRequiredParam<Real>("C_g", "C_g parameter for JWL eos");
  params.addRequiredParam<Real>("R_1g", "R_1g parameter for JWL eos");
  params.addRequiredParam<Real>("R_2g", "R_2g parameter for JWL eos");
  params.addRequiredParam<Real>("omega", "omega parameter for gas JWL eos");

  //parameters for JWL EOS for the solid phase
  
  params.addRequiredParam<Real>("A_s", "A_s parameter for solid JWL eos");
  params.addRequiredParam<Real>("B_s", "B_s parameter for solid JWL eos");
  params.addRequiredParam<Real>("C_s", "C_s parameter for solid JWL eos");
  params.addRequiredParam<Real>("R_1s", "R_1s parameter for solid JWL eos");
  params.addRequiredParam<Real>("R_2s", "R_2s parameter for solid JWL eos");
  params.addRequiredParam<Real>("omega_s", "omega parameter for solid JWL eos");

  //evolution equation parameters

  //params.addRequiredCoupledVar("lambda", "lambda"); //retrieve lambda as a coupled variable
  params.addRequiredParam<Real>("I_ev", "I parameter for evolution equation");
  params.addRequiredParam<Real>("G1", "G1 parameter for evolution equation");
  params.addRequiredParam<Real>("G2", "G2 parameter for evolution equation");
  params.addRequiredParam<Real>("x0", "x0 parameter for evolution equation");
  params.addRequiredParam<Real>("y0", "y0 parameter for evolution equation");
  params.addRequiredParam<Real>("x1", "x1 parameter for evolution equation");
  params.addRequiredParam<Real>("y1", "y1 parameter for evolution equation");
  params.addRequiredParam<Real>("z1", "z1 parameter for evolution equation");
  params.addRequiredParam<Real>("x2", "x2 parameter for evolution equation");
  params.addRequiredParam<Real>("y2", "y2 parameter for evolution equation");
  params.addRequiredParam<Real>("z2", "z2 parameter for evolution equation");
  params.addRequiredParam<Real>("a0", "a0 parameter for evolution equation");
  params.addRequiredParam<Real>("useArrhenius", "chemistry model to use");
  params.addRequiredParam<Real>("isTempDep", "is the JWL equation temperature dependent for the reacted products");
  params.addRequiredParam<Real>("isMG", "is the EOS for the solid phase MG or JWL");
  params.addRequiredParam<Real>("igmax", "maximum ignition phase ration value");
  params.addRequiredParam<Real>("gmax1", "maximum growth 1 value");
  return params;
}

ComputeMieGruneisenPlasticPFFractureStressDebugChem::ComputeMieGruneisenPlasticPFFractureStressDebugChem(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _Gamma(getParam<Real>("Gamma")),
    _density(getMaterialProperty<Real>("density")),
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _temperature(coupledValue("temperature")),
    _ref_temperature(getParam<Real>("reference_temperature")),
    _s(getParam<Real>("slope_UsUp")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _Le(getParam<Real>("Le")),
    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _pressure_eos(declareProperty<Real>("pressure_eos")),
    _plastic_strain(declareProperty<RankTwoTensor>(_base_name + "plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "plastic_strain")),
    _plastic_strain_rate(declareProperty<RankTwoTensor>(_base_name + "plastic_strain_rate")),
    _eqv_plastic_strain(declareProperty<Real>(_base_name + "eqv_plastic_strain")),
    _eqv_plastic_strain_rate(declareProperty<Real>(_base_name + "eqv_plastic_strain_rate")),
    _eqv_plastic_strain_rate_old(getMaterialPropertyOldByName<Real>(_base_name + "eqv_plastic_strain_rate")),
    _eqv_plastic_strain_old(getMaterialPropertyOld<Real>(_base_name + "eqv_plastic_strain")),
    _rotation_increment(getMaterialProperty<RankTwoTensor>(_base_name + "rotation_increment")),
    _rtol(getParam<Real>("rtol")),
    _ftol(getParam<Real>("ftol")),
    _eptol(getParam<Real>("eptol")),
    _deltaOuter(RankTwoTensor::Identity().times<i_, j_, k_, l_>(RankTwoTensor::Identity())),
    _deltaMixed(RankTwoTensor::Identity().times<i_, k_, j_, l_>(RankTwoTensor::Identity())),
    _W0p(declareProperty<Real>("W0p")),
    _stress_eos_elastic(declareProperty<RankTwoTensor>(_base_name + "stress_eos_elastic")),
    _stress_cpl_elastic(declareProperty<RankTwoTensor>(_base_name + "stress_cpl_elastic")),
    _stress_eos_plastic(declareProperty<RankTwoTensor>(_base_name + "stress_eos_plastic")),
    _stress_cpl_plastic(declareProperty<RankTwoTensor>(_base_name + "stress_cpl_plastic")),
    _A_yield(getParam<Real>("A_yield")),
    _B_yield(getParam<Real>("B_yield")),
    _C_yield(getParam<Real>("C_yield")),
    _a_corr_yield(getParam<Real>("a_corr_yield")),
    _m_yield(getParam<Real>("m_yield")),
    _n_yield(getParam<Real>("n_yield")),
    _strain_ref(getParam<Real>("strain_ref")),
    _T_melt_ref(getParam<Real>("T_melt_ref")),
    _maxiter(getParam<Real>("maxiter")),
    _toly(getParam<Real>("toly")),
    _k_alpha(getParam<Real>("k_alpha")),
    _thr_gas(getParam<Real>("thr_gas")),
    _yield_1(declareProperty<Real>(_base_name + "yield_1")),
    _yield_2(declareProperty<Real>(_base_name + "yield_2")),
    _yield_3(declareProperty<Real>(_base_name + "yield_3")),
    _yield_total(declareProperty<Real>(_base_name + "yield_total")),
    _dyield_dplastic(declareProperty<Real>(_base_name + "dyield_dplastic")),
    _theta_norm(declareProperty<Real>(_base_name + "theta_norm")),
    _T_melt(declareProperty<Real>(_base_name + "T_melt")),
    _V(declareProperty<Real>(_base_name + "V")),
    _thr_linear(getParam<Real>("thr_linear")),
    _eqv_trial(declareProperty<Real>(_base_name + "eqv_trial")),
    _yield3_thr(getParam<Real>(_base_name + "yield3_thr")),
    _dyield1(declareProperty<Real>(_base_name + "dyield1")),
    _ss_prop(declareProperty<Real>(_base_name + "ss_prop")),
    _gamma_pd(getParam<Real>("gamma_pd")),
    //chemical reaction parameters
    _chemY1(coupledValue("chemY1")),
    _chemY2(coupledValue("chemY2")),
    _chemY3(coupledValue("chemY3")),
    _chemY4(coupledValue("chemY4")),
    _chemY1_dot(declareProperty<Real>("chemY1_dot")),
    _chemY2_dot(declareProperty<Real>("chemY2_dot")),
    _chemY3_dot(declareProperty<Real>("chemY3_dot")),
    _chemY4_dot(declareProperty<Real>("chemY4_dot")),

    _Z1(getParam<Real>("Z1")),
    _Z2(getParam<Real>("Z2")),
    _Z3(getParam<Real>("Z3")),
    _Z4(getParam<Real>("Z4")),

    _E1(getParam<Real>("E1")),
    _E2(getParam<Real>("E2")),
    _E3(getParam<Real>("E3")),
    _R_const(getParam<Real>("R_const")),

    //gas phase parameters

    _A_g(getParam<Real>("A_g")),
    _B_g(getParam<Real>("B_g")),
    _C_g(getParam<Real>("C_g")),
    _R_1g(getParam<Real>("R_1g")),
    _R_2g(getParam<Real>("R_2g")),
    _omega(getParam<Real>("omega")),

    _A_s(getParam<Real>("A_s")),
    _B_s(getParam<Real>("B_s")),
    _C_s(getParam<Real>("C_s")),
    _R_1s(getParam<Real>("R_1s")),
    _R_2s(getParam<Real>("R_2s")),
    _omega_s(getParam<Real>("omega_s")),

    _p_g(declareProperty<Real>("p_g")),
    _p_s(declareProperty<Real>("p_s")),

    // chemical decomposition evolution equation parameters
    
    //_lambda_prop(declareProperty<Real>("lambda_prop")), //declares a material property to store the value of the lambda variable
    _lambda_dot(declareProperty<Real>("lambda_dot")),
    _lambda(declareProperty<Real>("lambda")),
    _lambda_old(getMaterialPropertyOld<Real>("lambda")),
    //_lambda_prop_old(getMaterialPropertyOld<Real>("lambda_prop")), //gets the previous value of lambda thru lambda_prop_old
    //_lambda(coupledValue("lambda")),
    _I_ev(getParam<Real>("I_ev")),
    _G1(getParam<Real>("G1")),
    _G2(getParam<Real>("G2")),

    _x0(getParam<Real>("x0")),
    _y0(getParam<Real>("y0")),
    _x1(getParam<Real>("x1")),
    _y1(getParam<Real>("y1")),
    _z1(getParam<Real>("z1")),
    _x2(getParam<Real>("x2")),
    _y2(getParam<Real>("y2")),
    _z2(getParam<Real>("z2")),
    _a0(getParam<Real>("a0")),

    _useArrhenius(getParam<Real>("useArrhenius")),
    _density_corr(declareProperty<Real>("density_corr")),
    _mu(declareProperty<Real>("mu")),

    _isTempDep(getParam<Real>("isTempDep")),
    _isMG(getParam<Real>("isMG")),
    _Cv_gas(getMaterialProperty<Real>("Cv_gas")), //defined in the material's parameters
    _Cv_solid(getMaterialProperty<Real>("Cv_solid")),
    _lambda1(declareProperty<Real>("lambda1")),
    _lambda2(declareProperty<Real>("lambda2")),
    _lambda3(declareProperty<Real>("lambda3")),
    _igmax(getParam<Real>("igmax")), //ignition evolution upper limit
    _gmax1(getParam<Real>("gmax1")), //growth1 evolution upper limit
    _pressure_eos_plastic(declareProperty<Real>("pressure_eos_plastic")),
    _jac_lambda(declareProperty<Real>("jac_lambda"))
{
}

void
ComputeMieGruneisenPlasticPFFractureStressDebugChem::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  _plastic_strain[_qp].zero();
  _eqv_plastic_strain[_qp] = 0.0;
  _eqv_plastic_strain_rate[_qp] = 0.0; //initialize the variable
  _plastic_strain_rate[_qp].zero();

  const double B=1E10;
  const double PI=3.141592653589793238463;

  //initialize the chemical composition

  //_lambda[_qp] = 0.0;
  //_lambda1[_qp] = 0.0;
  //_lambda2[_qp] = 0.0;
  //_lambda3[_qp] = 0.0;

 
}

void
ComputeMieGruneisenPlasticPFFractureStressDebugChem::computeQpStress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  //  Calculate "elastic" stresses due to equation of state
  // Calculate pressure from Mie Gruneisen (Menon, 2014), (Zhang, 2011)
  Real K0, delta, eta, temperature, p_s, ss_prop;
  RankTwoTensor stress_eos, stress, stress_cpl;
  K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  ss_prop = std::sqrt(K0 / _density[_qp]); //store sound speed as a function of actual bulk modulus from elasticity tensor
  _ss_prop[_qp] = ss_prop;
  _bulk_modulus[_qp] = K0;
  delta = _mechanical_strain[_qp].trace();
  eta = - delta;
  temperature = _temperature[_qp];
  Real V = 1.0 + _mechanical_strain[_qp].trace();
  _mu[_qp] = - delta / (1.0 + delta);
  _density_corr[_qp] = _density[_qp] * (1.0 / V);

  //eos for the unreacted phase (Y1): Mie-Gruneisen

  if (_isMG == 1) {
    p_s =  (- K0 * eta * (1.0 - (_Gamma * eta / 2.0)) / std::pow((1.0 - _s * eta), 2.0)) - (_Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature) * (1.0 / V));
  }
  else if (_isMG == 0) {
    p_s = (_A_s * std::exp(- _R_1s * (V))) + (_B_s * std::exp(- _R_2s * (V))) + (_omega_s * _Cv_solid[_qp] * _temperature[_qp]);// * (1.0 / V));
    p_s *= -1.0; //invert sign to account for compressive loading condition
  }
  
  _p_s[_qp] = p_s; //store on declared material property

  //eos for the reacted phases (1 - Y1 = Y2 + Y3 + Y4): Jones-Wilkins-Lee C_form

  if (_isTempDep == 1) {
    _p_g[_qp] = (_A_g * std::exp(- _R_1g * (V))) + (_B_g * std::exp(- _R_2g * (V))) + (_omega * _Cv_gas[_qp] * _temperature[_qp]);// * (1.0 / V));
  }
  else if (_isTempDep == 0) {
    _p_g[_qp] = (_A_g * std::exp(- _R_1g * (V))) + (_B_g * std::exp(- _R_2g * (V))) + (_C_g * std::pow((V), - (1.0 + _omega)));
  }

  _p_g[_qp] *= -1; //make gas pressure negative for keeping the same notation as the unreacted pressure values

  //combined phase eos  
  Real peos;
  if (_useArrhenius == 1) {
    peos = ((1.0 - _chemY1[_qp]) * _p_g[_qp]) + (_chemY1[_qp] * _p_s[_qp]);
  }
  else if (_useArrhenius == 0) {
    peos = ((1.0 - _lambda[_qp]) * _p_s[_qp]) + (_lambda[_qp] * _p_g[_qp]);
  }

  peos *= 1; //make negative for compression 

  //evolution equation for chemistry single step reaction
  //declare the conditions for each evolution stage
  //current rate of change of lambda

  Real frac = (1.0 - _lambda[_qp]);
  if (_lambda[_qp] <= _igmax) {
    _lambda1[_qp] = _I_ev * std::pow(frac, _x0) * std::pow(_mu[_qp] - _a0, _y0);
  }

  else {
    _lambda1[_qp] = 0.0;
  }

  if (_lambda[_qp] <= _gmax1) {
    _lambda2[_qp] = _G1 * std::pow(frac, _x1) * std::pow(_lambda[_qp], _y1) * std::pow(std::abs(peos), _z1);
  }

  else {
    _lambda2[_qp] = 0.0;
  }

  if (_lambda[_qp] >= _gmax1) {
    _lambda3[_qp] = _G2 * std::pow(frac, _x2) * std::pow(_lambda[_qp], _y2) * std::pow(std::abs(peos), _z2);
  }

  else {
    _lambda3[_qp] = 0.0;
  }

  _lambda_dot[_qp] = _lambda1[_qp] + _lambda2[_qp] + _lambda3[_qp];

  // compute the jacobian for the kernel formulation

  Real jac1, jac2, jac3;
  Real lambda;
  lambda = _lambda[_qp]; //internal variable

  jac1 = - _x0 * _I_ev * std::pow((1.0 - lambda), _x0 - 1.0) * std::pow(((1.0 / _V[_qp]) - 1.0 - _a0), _y0);

  jac2 = - _x1 * _G1 * std::pow((1.0 - lambda), _x1 - 1.0) * std::pow(lambda, _y1) * std::pow(peos, _z1);
  jac2 += _G1 * std::pow((1.0 - lambda), _x1) * _y1 * std::pow(lambda, _y1 - 1.0) * std::pow(peos, _z1);

  jac3 = - _x2 * _G2 * std::pow((1.0 - lambda), _x2 - 1.0) * std::pow(lambda, _y2) * std::pow(peos, _z2);
  jac3 += _G2 * std::pow((1.0 - lambda), _x2) * _y2 * std::pow(lambda, _y2 - 1.0) * std::pow(peos, _z2);

  _jac_lambda[_qp] = jac1 + jac2 + jac3; //compute the entire jacobian for the kernel formulation 

  // ------------------------
  //NOTE: this code computes the RHS of the evloution equation for lambda. It stores the RHS in a material
  //      property that is later retrieved by a kernel to compute the residual of the equation
  //      lambda_dot - f(lambda) = 0.0, where f(lambda) is the pressure-dependent evolution equation RHS.
  // ------------------------

   _lambda[_qp] = _lambda_old[_qp] + (_lambda_dot[_qp] * _dt);

  if (_lambda[_qp] >= 1.0) { //condition supplied to enforce lambda to be less than or equal to one

    _lambda[_qp] = 1.0; //fix the value of lambda at it's maximum possible value

  }

  //define the chemistry rates. stores the rate in material properties to then be used in the kernel

  _chemY1_dot[_qp] = - _chemY1[_qp] * _Z1 * std::exp(- _E1 / (_R_const * _temperature[_qp]));
  _chemY2_dot[_qp] = (_chemY1[_qp] * _Z1 * std::exp(- _E1 / (_R_const * _temperature[_qp]))) - (_chemY2[_qp] * _Z2 * std::exp(- _E2 / (_R_const * _temperature[_qp])));
  _chemY3_dot[_qp] = (_chemY2[_qp] * _Z2 * std::exp(- _E2 / (_R_const * _temperature[_qp]))) - ((_chemY3[_qp] * _chemY3[_qp]) * _Z3 * std::exp(- _E3 / (_R_const * _temperature[_qp])));
  _chemY4_dot[_qp] = (_chemY3[_qp] * _chemY3[_qp]) * _Z3 * std::exp(- _E3 / (_R_const * _temperature[_qp]));

  //update variables
  _pressure_eos[_qp] = peos;
  stress_eos = peos * I2;
  _stress_eos_elastic[_qp] = stress_eos;
  stress_cpl = _elasticity_tensor[_qp] * (_elastic_strain_old[_qp] + _strain_increment[_qp]) - K0 * delta * I2;
  _stress_cpl_elastic[_qp] = stress_cpl;
  stress = stress_eos + stress_cpl; //changed to account for the volumetric expansion stress
  _stress[_qp] = stress_eos; //preliminary stress value storing

  //  Calculate plastic strain update
  Real f, flow_incr, rep, dflow_incr, deqvpstrain, fq, yield_stress, err1, err2, err3, toly, eqvpstrain;
  RankTwoTensor flow_dirn, resid, ddsig, delta_dp, df_dsig, plastic_strain;
  RankFourTensor dr_dsig, dr_dsig_inv;

  flow_incr = 0.0;
  dflow_incr = 0.0;
  deqvpstrain = 0.0;
  toly = 1.0e-3;

  unsigned int iter = 0;
  unsigned int maxiter = 1000;

  // this variables are set to the previous timestep to then add the computed increment on the current timestep
  eqvpstrain = _eqv_plastic_strain_old[_qp];
  plastic_strain = _plastic_strain_old[_qp];
  _eqv_plastic_strain_rate[_qp] = _eqv_plastic_strain_rate_old[_qp]; //set the plastic strain rate to the previous value to update in current timestep CHECK

  //_eqv_plastic_strain_rate[_qp] = (_eqv_plastic_strain[_qp] - _eqv_plastic_strain_old[_qp]) / _dt;
  //_plastic_strain_rate[_qp] = (_plastic_strain[_qp] - _plastic_strain_old[_qp]) / _dt;

  yield_stress = getYieldStress(eqvpstrain, _eqv_plastic_strain_rate[_qp], _temperature[_qp]); // gets the yield stress at the old configuration but current temperature CHECK
  if (yieldFunction(stress, yield_stress) > _toly)
  {
    // the sig just calculated is inadmissable.  We must return to the yield surface.
    // This is done iteratively, using a Newton-Raphson process.
    delta_dp.zero();
    flow_dirn = flowPotential(stress); //flow direction at current EOS + CPL(deviatoric) stress config CHECK

    resid = flow_dirn * flow_incr - delta_dp; // Residual 1 - refer Hughes Simo
    f = yieldFunction(stress, yield_stress); // computes the yield difference (surface - stress) for the current stress and with the previous yield definition
    rep = -eqvpstrain + _eqv_plastic_strain_old[_qp] - flow_incr * internalPotential(); // Residual 3 rep=0

    err1 = resid.L2norm();
    err2 = std::abs(f);
    err3 = std::abs(rep);

    while ((err1 > _rtol || err2 > _ftol || err3 > _eptol) &&
           iter < _maxiter) // Stress update iteration (hardness fixed)
    {
      iter++;

      df_dsig = dyieldFunction_dstress(stress);
      getJac(stress, _elasticity_tensor[_qp], flow_incr, dr_dsig);   // gets dr_dsig = d(resid_ij)/d(sig_kl)
      fq = dyieldFunction_dinternal(eqvpstrain, _eqv_plastic_strain_rate[_qp], _temperature[_qp]); // d(f)/d(eqvpstrain) at the same configuration as the yield surface CHECK
      dr_dsig_inv = dr_dsig.invSymm();

      dflow_incr = (f - df_dsig.doubleContraction(dr_dsig_inv * resid) + fq * rep) / (df_dsig.doubleContraction(dr_dsig_inv * flow_dirn) - fq);
      ddsig = dr_dsig_inv * (-resid - flow_dirn * dflow_incr); // from solving the top row of linear system, given dflow_incr
      deqvpstrain = rep + dflow_incr; // from solving the bottom row of linear system, given dflow_incr

      // update the variables with the increment computed during the current timestep
      flow_incr += dflow_incr;
      delta_dp -= _elasticity_tensor[_qp].invSymm() * ddsig;
      stress += ddsig;
      eqvpstrain += deqvpstrain;

      // evaluate the RHS equations ready for next Newton-Raphson iteration
      flow_dirn = flowPotential(stress);
      resid = flow_dirn * flow_incr - delta_dp;
      f = yieldFunction(stress, yield_stress);
      rep = -eqvpstrain + _eqv_plastic_strain_old[_qp] - flow_incr * internalPotential();

      err1 = resid.L2norm();
      err2 = std::abs(f);
      err3 = std::abs(rep);
    }
    

    if (iter >= _maxiter)
      mooseError("Constitutive failure");

    plastic_strain += delta_dp;
    
  }
  //store the variables for visualization purposes

  _eqv_plastic_strain[_qp] = eqvpstrain; //this already has the increment added
  // _eqv_plastic_strain_rate[_qp] = deqvpstrain / _dt; // add to the old strain rate the increment over the timestep CHECK
  _eqv_plastic_strain_rate[_qp] = (_eqv_plastic_strain[_qp] - _eqv_plastic_strain_old[_qp]) / _dt; // add to the old strain rate the increment over the timestep CHECK
  
  // Rotate plastic strain tensor to the current configuration
  _plastic_strain[_qp] = plastic_strain;
  _plastic_strain_rate[_qp] = (_plastic_strain[_qp] - _plastic_strain_old[_qp]) / _dt;
  
  RankTwoTensor p_rate = _plastic_strain_rate[_qp];
  unsigned int i, j, k; //define the orthogonal basis vector for iterating over the entries of the tensor
  for (i = 0; i < 3; i++) {
  	for (j = 0; j < 3; j++) {
    		if (p_rate(i, j) > _gamma_pd) {
    			p_rate(i, j) = _gamma_pd;
    		}
    		else if (p_rate(i, j) < - _gamma_pd) {
    			p_rate(i, j) = - _gamma_pd;
    		}
    	}
  }
  _plastic_strain_rate[_qp] = p_rate; //update the value of the tensor with the iterated variable
  
  //limit the equivalent plastic strain rate
  
  if (_eqv_plastic_strain_rate[_qp] >= _gamma_pd) {
  	_eqv_plastic_strain_rate[_qp] = _gamma_pd;
  }

  _eqv_trial[_qp] = std::sqrt((2.0 /3.0) * (_plastic_strain[_qp].doubleContraction(_plastic_strain[_qp]))); // to compare algorithm's calculation with analitical definition

  // Calculate the elastic strain_increment
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _plastic_strain[_qp];

  //  Store plasticity updated stress
  peos = stress.trace()/3;
  _pressure_eos_plastic[_qp] = peos;
  stress_eos = peos * I2;
  _stress_eos_plastic[_qp] = stress_eos;
  stress_cpl = stress - stress_eos;
  _stress_cpl_plastic[_qp] = stress_cpl;
  _stress[_qp] = stress; //CHECK the additive formulation of this stress decomposition

  // Calculate bulk-viscosity stress term
  Real trD, jacob, q_bv;
  trD = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  jacob = 1.0 + _mechanical_strain[_qp].trace();
  q_bv = 0.0;
  if (jacob < 1.0) {
    q_bv = ( _C0 * _density[_qp] * trD * std::abs(trD) * std::pow(_Le,2.0) / std::pow(jacob,2.0) ) + ( _C1 * _density[_qp] * ss_prop * trD * _Le / jacob ); //edited sound speed CHECK
  }
  _stress[_qp] += q_bv * I2;

  // //chemical reaction evolution functions
  // if (_useArrhenius == 1) {
  //   _chemY1_dot[_qp] = chem_react_1(_chemY1[_qp], _temperature[_qp]);
  //   _chemY2_dot[_qp] = chem_react_2(_chemY1[_qp], _chemY1_dot[_qp], _chemY2[_qp], _temperature[_qp]);
  //   _chemY3_dot[_qp] = chem_react_3(_chemY2[_qp], _chemY2_dot[_qp], _chemY3[_qp], _temperature[_qp]);
  //   _chemY4_dot[_qp] = chem_react_4(_chemY3[_qp], _chemY3_dot[_qp], _temperature[_qp]);

  //   //formulate the explicit integration to obtain each species

  //   _chemY1[_qp] = _chemY1_old[_qp] + _chemY1_dot[_qp] * _dt;
  //   _chemY2[_qp] = _chemY2_old[_qp] + _chemY2_dot[_qp] * _dt;
  //   _chemY3[_qp] = _chemY3_old[_qp] + _chemY3_dot[_qp] * _dt;
  //   _chemY4[_qp] = _chemY4_old[_qp] + _chemY4_dot[_qp] * _dt;
  // }
}

Real
ComputeMieGruneisenPlasticPFFractureStressDebugChem::yieldFunction(const RankTwoTensor & stress, const Real yield_stress)
{
  return getSigEqv(stress) - yield_stress;
}

RankTwoTensor
ComputeMieGruneisenPlasticPFFractureStressDebugChem::dyieldFunction_dstress(const RankTwoTensor & sig)
{
  RankTwoTensor deriv = sig.dsecondInvariant();
  deriv *= std::sqrt(3.0 / sig.secondInvariant()) / 2.0;
  return deriv;
}

Real
ComputeMieGruneisenPlasticPFFractureStressDebugChem::dyieldFunction_dinternal(const Real eqvpstrain, const Real eqv_plastic_strain_rate, const Real Temperature)
{
  return -getdYieldStressdPlasticStrain(eqvpstrain, eqv_plastic_strain_rate, Temperature);
}

RankTwoTensor
ComputeMieGruneisenPlasticPFFractureStressDebugChem::flowPotential(const RankTwoTensor & sig)
{
  return dyieldFunction_dstress(sig); // this plasticity model assumes associative flow
}

Real
ComputeMieGruneisenPlasticPFFractureStressDebugChem::internalPotential()
{
  return -1;
}

Real
ComputeMieGruneisenPlasticPFFractureStressDebugChem::getSigEqv(const RankTwoTensor & stress)
{
  return std::sqrt(3 * stress.secondInvariant());
}

// Jacobian for stress update algorithm
void
ComputeMieGruneisenPlasticPFFractureStressDebugChem::getJac(const RankTwoTensor & sig,
                                    const RankFourTensor & E_ijkl,
                                    Real flow_incr,
                                    RankFourTensor & dresid_dsig)
{
  RankTwoTensor sig_dev, df_dsig, flow_dirn;
  RankTwoTensor dfi_dft, dfi_dsig;
  RankFourTensor dft_dsig, dfd_dft, dfd_dsig;
  Real sig_eqv;
  Real f1, f2, f3;
  RankFourTensor temp;

  sig_dev = sig.deviatoric();
  sig_eqv = getSigEqv(sig);
  df_dsig = dyieldFunction_dstress(sig);
  flow_dirn = flowPotential(sig);

  f1 = 3.0 / (2.0 * sig_eqv);
  f2 = f1 / 3.0;
  f3 = 9.0 / (4.0 * Utility::pow<3>(sig_eqv));

  dft_dsig = f1 * _deltaMixed - f2 * _deltaOuter - f3 * sig_dev.outerProduct(sig_dev);

  dfd_dsig = dft_dsig;
  dresid_dsig = E_ijkl.invSymm() + dfd_dsig * flow_incr;
}

// Obtain yield stress for a given equivalent plastic strain (input)
Real
ComputeMieGruneisenPlasticPFFractureStressDebugChem::getYieldStress(const Real eqvpstrain, const Real eqv_plastic_strain_rate, const Real Temperature)
{
  Real sigma_Y_1, sigma_Y_2, sigma_Y_3, sigma_Y, theta, T_melt;
  // _eqv_plastic_strain_rate[_qp] = (_eqv_plastic_strain[_qp] - _eqv_plastic_strain_old[_qp]) / _dt;

  sigma_Y = 0.0;
  sigma_Y_1 = _A_yield + (_B_yield * std::pow(eqvpstrain, _n_yield));
  _yield_1[_qp] = sigma_Y_1; //store on material property
  sigma_Y += sigma_Y_1;
  
  //sigma_Y_2 = 1.0 + (_C_yield * std::max(0.0, std::log(eqv_plastic_strain_rate / _strain_ref)));
  
  if (eqv_plastic_strain_rate < _strain_ref) {
    sigma_Y_2 = 1.0;
  }
  else {
    sigma_Y_2 = 1.0 + (_C_yield * std::log(eqv_plastic_strain_rate / _strain_ref));
  }
  
  sigma_Y *= sigma_Y_2;
  _yield_2[_qp] = sigma_Y_2;

  // compute temperature related terms
  Real V, exp_t; // where V = v/v_0
  exp_t = (2.0 / 3.0) * (_Gamma - _a_corr_yield - 1.0);
  V = 1 + _mechanical_strain[_qp].trace();
  _V[_qp] = V; //store relative volume change on declared material variable
  
  T_melt = _T_melt_ref * std::exp(2.0 * _a_corr_yield * (1.0 - V)) * std::pow((1.0 / V), exp_t);
  theta = (Temperature - _ref_temperature) / (T_melt - _ref_temperature);
  
  if (theta >= _yield3_thr) {
  	theta = _yield3_thr;
  }
    
  _T_melt[_qp] = T_melt; //store on material property
  _theta_norm[_qp] = theta; //store on material property

  sigma_Y_3 = (1.0 - std::pow(theta, _m_yield));
  
  _yield_3[_qp] = sigma_Y_3;
  sigma_Y *= sigma_Y_3; // final form of yield stress
  _yield_total[_qp] = sigma_Y; //store whole expression yield stress
  return sigma_Y;
}

Real
ComputeMieGruneisenPlasticPFFractureStressDebugChem::getdYieldStressdPlasticStrain(const Real eqvpstrain, const Real eqv_plastic_strain_rate, const Real Temperature)
{
  Real sigma_Y_1, sigma_Y_2, sigma_Y_3, sigma_Y, theta, T_melt;
  Real dy_dplastic, dy_dplastic_1, dy_dplastic_2;
  // _eqv_plastic_strain_rate[_qp] = (_eqv_plastic_strain[_qp] - _eqv_plastic_strain_old[_qp]) / _dt;
  
  dy_dplastic = 0.0;
  dy_dplastic_1 = 0.0;
  _dyield1[_qp] = 0.0;
  //_dyield_dplastic[_qp] = 0.0;

  Real V, exp_t, thr_linear;
  thr_linear = 0.25; //threshold for avoiding singularity on the yield stress derivative definition
  exp_t = (2.0 / 3.0) * ((3.0 * _Gamma) - (3.0 * _a_corr_yield) - 1.0); //this is a constant value
  V = 1 + _mechanical_strain[_qp].trace(); // as per definition of eta = 1 - v/v0

  T_melt = _T_melt_ref * std::exp(2.0 * _a_corr_yield * (1 - V)) * std::pow((1.0 / V), exp_t);
  theta = (Temperature - _ref_temperature) / (T_melt - _ref_temperature);
  
  // introduce threshold to ensure thermal softening doesn't affect convergence
  if (theta >= _yield3_thr) {
  	theta = _yield3_thr;
  }
  
  //fixing singularity issues
  if (eqv_plastic_strain_rate < _strain_ref) {
      if (eqvpstrain < _thr_linear) {
          dy_dplastic_1 += _B_yield * _n_yield * std::pow(_thr_linear, (_n_yield - 1.0))
                     * 1.0
                     * (1.0 - std::pow(theta, _m_yield));
          _dyield1[_qp] += _B_yield * _n_yield * std::pow(_thr_linear, (_n_yield - 1.0));
      }
      else {
          dy_dplastic_1 += _B_yield * _n_yield * std::pow(eqvpstrain, (_n_yield - 1.0))
                     * 1.0
                     * (1.0 - std::pow(theta, _m_yield));
          _dyield1[_qp] += _B_yield * _n_yield * std::pow(eqvpstrain, (_n_yield - 1.0));          
      }
  }
  else {
      if (eqvpstrain < _thr_linear) {
          dy_dplastic_1 += _B_yield * _n_yield * std::pow(_thr_linear, (_n_yield - 1.0))
                     * (1.0 + _C_yield * std::log(eqv_plastic_strain_rate / _strain_ref))
                     * (1.0 - std::pow(theta, _m_yield));
      }
      else {
          dy_dplastic_1 += _B_yield * _n_yield * std::pow(eqvpstrain, (_n_yield - 1.0))
                     * (1.0 + _C_yield * std::log(eqv_plastic_strain_rate / _strain_ref))
                     * (1.0 - std::pow(theta, _m_yield));          
      }
  }
  
  //dy_dplastic_2 = (_A_yield + _B_yield * std::pow(std::abs(eqvpstrain), _n_yield))
  //                   * (_C_yield / eqv_plastic_strain_rate);
  //                   * (1.0 - std::pow(theta, _m_yield));
  
  dy_dplastic = dy_dplastic_1; // + dy_dplastic_2;
  _dyield_dplastic[_qp] += dy_dplastic; //store on material property
  return dy_dplastic;
                     
}

// Real
// ComputeMieGruneisenPlasticPFFractureStressDebugChem::chem_react_1(const Real Y1, const Real T) {
//   Real Y1_rate;
//   Y1_rate = - Y1 * _Z1 * std::exp(- (_E1) / (_R_const * T));
//   return Y1_rate;
// }

// Real
// ComputeMieGruneisenPlasticPFFractureStressDebugChem::chem_react_2(const Real Y1, const Real Y1_rate, const Real Y2, const Real T) {
//   Real Y2_rate;
//   Y2_rate = (Y1 * _Z1 * std::exp(- (_E1) / (_R_const * T))) - (Y2 * _Z2 * std::exp(- (_E2) / (_R_const * T)));
//   return Y2_rate;
// }

// Real 
// ComputeMieGruneisenPlasticPFFractureStressDebugChem::chem_react_3(const Real Y2, const Real Y2_rate, const Real Y3, const Real T) {
//   Real Y3_rate;
//   Y3_rate = (Y2 * _Z2 * std::exp(- (_E2) / (_R_const * T))) - ((Y3 * Y3) * _Z3 * std::exp(- (_E3) / (_R_const * T)));
//   return Y3_rate;
// }

//Real 
//ComputeMieGruneisenPlasticPFFractureStressDebugChem::chem_react_4(const Real Y3, const Real Y3_rate, const Real T) {
//  Real Y4_rate;
//  Y4_rate = (Y3 * Y3) * _Z3 * std::exp(- (_E3) / (_R_const * T));
//  return Y4_rate;
//}