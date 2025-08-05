/// Computes EoS stress with artificial viscosity and J2 plasticity with Johnson - Cook model for yield. Radial Return algorithm for returning to the yield surface
/// SGZ 2024

#include "ComputeMieGruneisenPlasticPFFractureStressDebug.h"

registerMooseObject("TensorMechanicsApp", ComputeMieGruneisenPlasticPFFractureStressDebug);

InputParameters
ComputeMieGruneisenPlasticPFFractureStressDebug::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addRequiredCoupledVar("c", "Name of damage variable");
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field fracture model, with small strain"
                             "Considers Mie Gruneisen EOS and artificial viscosity damping");
  params.addParam<bool>("use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<bool>("use_snes_vi_solver",false,"Use PETSc's SNES variational inequalities solver to enforce damage "
                        "irreversibility condition and restrict damage value <= 1.");
  params.addParam<MaterialPropertyName>("barrier_energy", "Name of material property for fracture energy barrier.");
  params.addParam<MaterialPropertyName>("E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>("D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<MaterialPropertyName>("I_name", "indicator", "Name of material property for damage indicator function.");
  params.addParam<MaterialPropertyName>("F_name", "local_fracture_energy", "Name of material property for local fracture energy function.");
  params.addRequiredParam<Real>("Gamma", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addParam<MaterialPropertyName>("density", "density", "Name of Material Property that provides the density");
  params.addParam<MaterialPropertyName>("specific_heat", "specific_heat", "Name of Material Property that provides the density");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("slope_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("Le","Maximum element size");
  // params.addRequiredParam<Real>("sound_speed","Speed of sound in the material");
  // params.addRequiredParam<std::vector<Real>>("yield_stress","Input data as pairs of equivalent plastic strain and yield stress: Should start with equivalent plastic strain 0");
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
  return params;
}

ComputeMieGruneisenPlasticPFFractureStressDebug::ComputeMieGruneisenPlasticPFFractureStressDebug(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _c(coupledValue("c")),
    _l(getMaterialProperty<Real>("l")),
    _pressure(getDefaultMaterialProperty<Real>("fracture_pressure")),
    _gc(getMaterialProperty<Real>("gc_prop")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _use_snes_vi_solver(getParam<bool>("use_snes_vi_solver")),
    _H(declareProperty<Real>("hist")),
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _barrier(getDefaultMaterialProperty<Real>("barrier_energy")),
    _E(declareProperty<Real>(getParam<MaterialPropertyName>("E_name"))),
    _dEdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name())),
    _d2Ed2c(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _D(getMaterialProperty<Real>("D_name")),
    _dDdc(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name())),
    _d2Dd2c(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _I(getDefaultMaterialProperty<Real>("I_name")),
    _dIdc(getMaterialPropertyDerivative<Real>("I_name", getVar("c", 0)->name())),
    _d2Id2c(getMaterialPropertyDerivative<Real>("I_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
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
    // _sound_speed(getParam<Real>("sound_speed")),
    _bulk_modulus(declareProperty<Real>("bulk_modulus")),
    _pressure_eos(declareProperty<Real>("pressure_eos")),
    // _yield_stress_vector(getParam<std::vector<Real>>("yield_stress")), // Read from input file
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
    _gamma_pd(getParam<Real>("gamma_pd"))
{
}

void
ComputeMieGruneisenPlasticPFFractureStressDebug::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  _plastic_strain[_qp].zero();
  _eqv_plastic_strain[_qp] = 0.0;
  _eqv_plastic_strain_rate[_qp] = 0.0; //initialize the variable
  _plastic_strain_rate[_qp].zero();

  _H[_qp] = 0.0;
  const double B=1E10;
  const double PI=3.141592653589793238463;

 //if(std::abs(_q_point[_qp](0)-12.5)<2){if(std::abs(_q_point[_qp](1)-0)<(0.3*_l[_qp])){_H[_qp] = B*(_gc[_qp]/4/(0.5*_l[_qp]))* std::exp(-std::abs(_q_point[_qp](0)-12.5)/(0.3*_l[_qp])); } }
}

void
ComputeMieGruneisenPlasticPFFractureStressDebug::computeQpStress()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  //  Calculate "elastic" stresses due to equation of state
  // Calculate pressure from Mie Gruneisen (Menon, 2014), (Zhang, 2011)
  // https://en.wikipedia.org/wiki/Mie%E2%80%93Gruneisen_equation_of_state
  Real K0, delta, eta, temperature, peos, ss_prop;
  RankTwoTensor stress_eos, stress, stress_cpl;
  K0 = (1.0 / 9.0) * I2.doubleContraction(_elasticity_tensor[_qp] * I2);
  ss_prop = std::sqrt(K0 / _density[_qp]); //store sound speed as a function of actual bulk modulus from elasticity tensor
  _ss_prop[_qp] = ss_prop;
  _bulk_modulus[_qp] = K0;
  delta = _mechanical_strain[_qp].trace();
  eta = - delta;
  temperature = _temperature[_qp];
  peos = - K0 * eta * (1.0 - (_Gamma * eta / 2.0)) / std::pow((1.0 - _s * eta), 2.0) - _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature) * eta; //correct formulation according to MieGruneisen EOS
  
  if (_c[_qp] > _thr_gas) {
    peos += _c[_qp] * _k_alpha * _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature) * eta; //gas expansion pressure term dissabled
  }
  
  _pressure_eos[_qp] = peos;
  stress_eos = peos * I2;
  _stress_eos_elastic[_qp] = stress_eos;
  stress_cpl = _elasticity_tensor[_qp] * (_elastic_strain_old[_qp] + _strain_increment[_qp]) - K0 * delta * I2;
  _stress_cpl_elastic[_qp] = stress_cpl;
  stress = stress_eos + stress_cpl;

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
  unsigned int i, j, k;
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
  //set the phonon drag limit for the shear strain rate
  
  //RankTwoTensor gamma_pd_tensor = _gamma_pd * std::ones(2);
  
  //RankTwoTensor p_rate = _plastic_strain_rate[_qp];
  
  //unsigned int i, j;
  //for (i = 0, j = 0; i < 3, j < 3; i++, j++){
  //  if (p_rate(i, j) >= _gamma_pd) {
  //    p_rate(i, j) = _gamma_pd;
  //    _plastic_strain_rate[_qp] = p_rate;
  //  }
  //}
  
  //limit the equivalent plastic strain rate
  
  if (_eqv_plastic_strain_rate[_qp] >= _gamma_pd) {
  	_eqv_plastic_strain_rate[_qp] = _gamma_pd;
  }

  
  _eqv_trial[_qp] = std::sqrt((2.0 /3.0) * (_plastic_strain[_qp].doubleContraction(_plastic_strain[_qp]))); // to compare algorithm's calculation with analitical definition

  // Calculate the elastic strain_increment
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _plastic_strain[_qp];

  //  Store plasticity updated stress
  peos = stress.trace()/3;
  stress_eos = peos * I2;
  _stress_eos_plastic[_qp] = stress_eos;
  stress_cpl = stress - stress_eos;
  _stress_cpl_plastic[_qp] = stress_cpl;

  // Create the positive and negative projection tensors
  std::vector<Real> eigval;
  RankTwoTensor eigvec;
  RankFourTensor Ppos = stress_cpl.positiveProjectionEigenDecomposition(eigval, eigvec);

  // Project the positive and negative stresses
  Real peos_pos;
  RankTwoTensor stress_eos_pos, stress_eos_neg, stress_cpl_pos, stress_cpl_neg, stress0pos, stress0neg;
  peos_pos = (std::abs(peos) + peos) / 2.0;
  stress_eos_pos = peos_pos * I2;
  stress_eos_neg = stress_eos - stress_eos_pos;
  stress_cpl_pos = Ppos * stress_cpl;
  stress_cpl_neg = stress_cpl - stress_cpl_pos;
  stress0pos = stress_eos_pos + stress_cpl_pos;
  stress0neg = stress - stress0pos;

  // Compute the positive and negative elastic energies
  Real F_pos, F_neg;
  Real A, B, C;
  A = K0 * (_Gamma * (1.0/(2.0 * std::pow(_s,2.0)) + 0.5 - 1.0 / _s) - 1.0);
  B = K0 * (_Gamma * (1.0 / _s - 0.5) + 1.0);
  C = - K0 * _Gamma / (2.0 * std::pow(_s,2.0));
  delta = _mechanical_strain[_qp].trace();
  eta = - _elastic_strain[_qp].trace();
  temperature = _temperature[_qp];
  if (delta>=0.0) {
    F_pos = A / (_s - std::pow(_s,2.0) * eta) - B * std::log(1.0 - _s * eta) / _s + C * eta - A / _s -
            _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature) * eta;
    F_neg = 0.0;
  }
  else {
    F_pos = 0.0;
    F_neg = A / (_s - std::pow(_s,2.0) * eta) - B * std::log(1.0 - _s * eta) / _s + C * eta - A / _s -
            _Gamma * _density[_qp] * _specific_heat[_qp] * (temperature - _ref_temperature) * eta;
  }

  F_pos += (stress_cpl_pos).doubleContraction(_elastic_strain[_qp]) / 2.0;
  F_neg += (stress_cpl_neg).doubleContraction(_elastic_strain[_qp]) / 2.0;

  _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];

  _Jacobian_mult[_qp] = (I4sym - (1 - _D[_qp]) * Ppos) * _elasticity_tensor[_qp];

  // Calculate bulk-viscosity stress term
  Real trD, jacob, q_bv;
  trD = (_mechanical_strain[_qp].trace() - _mechanical_strain_old[_qp].trace()) / _dt;
  jacob = 1.0 + _mechanical_strain[_qp].trace();
  q_bv = 0.0;
  if (jacob < 1.0) {
    q_bv = ( _C0 * _density[_qp] * trD * std::abs(trD) * std::pow(_Le,2.0) / std::pow(jacob,2.0) ) + ( _C1 * _density[_qp] * ss_prop * trD * _Le / jacob ); //edited sound speed CHECK
  }
  _stress[_qp] += q_bv * I2;


  // // Assign history variable
  Real hist_variable = _H_old[_qp];
  if (_use_snes_vi_solver)
  {
    _H[_qp] = F_pos;

    if (_use_current_hist)
      hist_variable = _H[_qp];
  }
  else
  {
    if (F_pos > _H_old[_qp])
      _H[_qp] = F_pos;
    else
      _H[_qp] = _H_old[_qp];

    if (_use_current_hist)
      hist_variable = _H[_qp];

    if (hist_variable < _barrier[_qp])
      hist_variable = _barrier[_qp];
  }

  // Elastic free energy density
  _E[_qp] = hist_variable * _D[_qp] + F_neg - _pressure[_qp] * _elastic_strain[_qp].trace() * _I[_qp];
  _dEdc[_qp] = hist_variable * _dDdc[_qp] - _pressure[_qp] * _elastic_strain[_qp].trace() * _dIdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp] - _pressure[_qp] * _elastic_strain[_qp].trace() * _d2Id2c[_qp];
}

Real
ComputeMieGruneisenPlasticPFFractureStressDebug::yieldFunction(const RankTwoTensor & stress, const Real yield_stress)
{
  return getSigEqv(stress) - yield_stress;
}

RankTwoTensor
ComputeMieGruneisenPlasticPFFractureStressDebug::dyieldFunction_dstress(const RankTwoTensor & sig)
{
  RankTwoTensor deriv = sig.dsecondInvariant();
  deriv *= std::sqrt(3.0 / sig.secondInvariant()) / 2.0;
  return deriv;
}

Real
ComputeMieGruneisenPlasticPFFractureStressDebug::dyieldFunction_dinternal(const Real eqvpstrain, const Real eqv_plastic_strain_rate, const Real Temperature)
{
  return -getdYieldStressdPlasticStrain(eqvpstrain, eqv_plastic_strain_rate, Temperature);
}

RankTwoTensor
ComputeMieGruneisenPlasticPFFractureStressDebug::flowPotential(const RankTwoTensor & sig)
{
  return dyieldFunction_dstress(sig); // this plasticity model assumes associative flow
}

Real
ComputeMieGruneisenPlasticPFFractureStressDebug::internalPotential()
{
  return -1;
}

Real
ComputeMieGruneisenPlasticPFFractureStressDebug::getSigEqv(const RankTwoTensor & stress)
{
  return std::sqrt(3 * stress.secondInvariant());
}

// Jacobian for stress update algorithm
void
ComputeMieGruneisenPlasticPFFractureStressDebug::getJac(const RankTwoTensor & sig,
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
ComputeMieGruneisenPlasticPFFractureStressDebug::getYieldStress(const Real eqvpstrain, const Real eqv_plastic_strain_rate, const Real Temperature)
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
ComputeMieGruneisenPlasticPFFractureStressDebug::getdYieldStressdPlasticStrain(const Real eqvpstrain, const Real eqv_plastic_strain_rate, const Real Temperature)
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
  exp_t = (2.0 / 3.0) * (_Gamma - _a_corr_yield - 1.0); //this is a constant value
  V = 1 + _mechanical_strain[_qp].trace(); // as per definition of eta = 1 - v/v0

  T_melt = _T_melt_ref * std::exp(2.0 * _a_corr_yield * (1 - V)) * std::pow((1.0 / V), exp_t);
  theta = (Temperature - _ref_temperature) / (T_melt - _ref_temperature);
  
  // introduce threshold to ensure thermal softening doesn't affect convergence
  if (theta >= _yield3_thr) {
  	theta = _yield3_thr;
  }
  
  //dy_dplastic_1 = _B_yield * _n_yield * std::pow(std::abs(eqvpstrain), (_n_yield - 1.0));
  
  //dy_dplastic_1 = _B_yield * _n_yield * std::pow((eqv_plastic_strain), (_n_yield - 1.0))
  //	         * (1.0 + _C_yield * std::max(0.0, std::log(_eqv_plastic_strain_rate[_qp] / _strain_ref)))
  //	         * (1.0 - std::pow(theta, _m_yield));
  
  
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
