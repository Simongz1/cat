E = 25.14
nu = 0.24

[GlobalParams]
  large_kinematics = true
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./dirac_switch_shock]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dirac_switch_react]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [disp_x]
    order = FIRST
  []
  [disp_y]
    order = FIRST
  []
  [temperature]
    order = FIRST
  []

	[./Y1]
    order = FIRST
    family = LAGRANGE
	[../]
	[./Y2]
    order = FIRST
    family = LAGRANGE
	[../]
	[./Y3]
    order = FIRST
    family = LAGRANGE
	[../]
[]

[AuxVariables]
	##dynamics variables
	[./vx]
    order = FIRST
    family = LAGRANGE
	[../]
	[./ax]
    order = FIRST
    family = LAGRANGE
	[../]
	[./vy]
    order = FIRST
    family = LAGRANGE
	[../]
	[./ay]
    order = FIRST
    family = LAGRANGE
	[../]

[]

[Mesh]
	type = GeneratedMesh
	dim = 2
	nx = 250
	ny = 30
	xmax = 475 #1.9mm
	xmin = 0
	ymax = 33 #1.1mm
	ymin = 0
	boundary_id = '0 1 2 3'
	boundary_name = 'bottom right top left'
  elem_type = QUAD4
[]

[Kernels]

  [dirac_switch_shock_time_derivative]
    type = ADTimeDerivative
    variable = dirac_switch_shock
  []
  [dirac_switch_react_time_derivative]
    type = ADTimeDerivative
    variable = dirac_switch_react
  []
  [rate_shock]
    type = MaskedBodyForce
    variable = dirac_switch_shock
    mask = rate_dirac_shock
  []
  [rate_react]
    type = MaskedBodyForce
    variable = dirac_switch_react
    mask = rate_dirac_react
  []
	
  ##temperature kernels
  [dTdt]
  	type = MassLumpedTimeDerivative
    variable = temperature
  []
  [nabla2T]
  	type = ADHeatConduction
    variable = temperature
    thermal_conductivity = 'k_corrected'
  []
  [sdx]
    type = TotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y'
  []
  [sdy]
    type = TotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y'
  []

	[./inertia_x]
		type = ADInertialForce
		variable = disp_x
		velocity = vx
		acceleration = ax
		beta = 0.3025 ###from dandekar 2019 sec 2.1
		gamma = 0.6 ###from dandekar 2019 sec 2.1
	[../]
	[./inertia_y]
		type = ADInertialForce
		variable = disp_y
		velocity = vy
		acceleration = ay
		beta = 0.3025 ###from dandekar 2019 sec 2.1
		gamma = 0.6 ###from dandekar 2019 sec 2.1
	[../]

  [chem_from_prop]
    type = ADMatHeatSource
    variable = temperature
    material_property = q_decomposition
    #displacements = 'disp_x disp_y'
    scalar = 1 #can switch sign
  []
  [PlasticHeat]
    type = ADMatHeatSource
    variable = temperature
    material_property = q_plastic
    #displacements = 'disp_x disp_y'
    scalar = 1 #to assign sign
  []
  [PressureHeat]
    type = ADMatHeatSource
    variable = temperature
    material_property = q_elastic
    #displacements = 'disp_x disp_y'
    scalar = 1
  []
  

  [shock_heat]
  	type = ADMISTERnetHeatShockLUMP
    variable = temperature
  []
  [shock_react]
  	type = ADMISTERnetHeatReactLUMP
    variable = temperature
  []

	[./dY1_dt]
		type = ADTimeDerivative #WORKS WITH ADTimeDerivative
		variable = Y1
	[../]
	[./dY2_dt]
		type = ADTimeDerivative
		variable = Y2
	[../]
	[./dY3_dt]
		type = ADTimeDerivative
		variable = Y3
	[../]

	[./Y1dot]
		type = ADY1_dot_RDX
		variable = Y1
	[../]
	[./Y2dot]
		type = ADY2_dot_RDX
		variable = Y2
		Y1 = Y1
	[../]
	[./Y3dot]
		type = ADY3_dot_RDX
		variable = Y3
		Y2 = Y2
	[../]

  #Test: surrogate chemistry sources
  [Y1_surr]
    type = ADSurrogateY
    variable = Y1
    surrogate_rate_name = Y3_dot_surrogate
    positive = true
  []
  [Y2_surr]
    type = ADSurrogateY
    variable = Y2
    surrogate_rate_name = Y2_dot_surrogate
    positive = true
  []
  [Y3_surr]
    type = ADSurrogateY
    variable = Y3
    surrogate_rate_name = Y3_dot_surrogate
    positive = false
  []
[]

[AuxKernels]
	[./vx]
	    	type = NewmarkVelAux
	    	variable = vx
	    	acceleration = ax
	    	gamma = 0.6
	[../]
	[./vy]
	    	type = NewmarkVelAux
	    	variable = vy
	    	acceleration = ay
	    	gamma = 0.6
	[../]

	[./ax]
	    	type = NewmarkAccelAux
	    	variable = ax
	    	displacement = disp_x
	    	velocity = vx
	    	beta = 0.3025
	[../]
	[./ay]
	    	type = NewmarkAccelAux
	    	variable = ay
	    	displacement = disp_y
	    	velocity = vy
	    	beta = 0.3025
	[../]	
[]

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  []
  [fix_y]
  	type = DirichletBC
  	variable = disp_y
  	boundary = 'top bottom'
  	value = 0.0
  []
  [./velocity]
  	type = PresetVelocity
  	boundary = right
  	velocity = -1.7
  	variable = disp_x
  [../]
[]

[ICs]
  [./dirac_switch_shock_ic]
    type = ConstantIC
    variable = dirac_switch_shock
    value = 0.0
  [../]
  [./dirac_switch_react_ic]
    type = ConstantIC
    variable = dirac_switch_react
    value = 0.0
  [../]
	[temp_IC]
		type = ConstantIC
		variable = temperature
		value = 300
	[]
	[./chemY1_IC]
		type = ConstantIC
		variable = Y1
		value = 1.0
	[../]
	[./chemY2_IC]
		type = ConstantIC
		variable = Y2
		value = 0.0
	[../]
	[./chemY3_IC]
		type = ConstantIC
		variable = Y3
		value = 0.0
	[../]
  [./RandomID]
    variable = density_i
    type = RandomIC
    min = 0
    max = 92
    seed = 1000412001
  [../]
  #[porosity]
  #  type = MultiSmoothCircleIC
  #  variable = density_i
  #  invalue = 90
  #  outvalue = 1
  #  bubspac = 15
  #  numbub = 10
  #  radius = 0.95
  #[]
[]

[Functions]
  [./timestep]
    type = ParsedFunction
    value = 'if(t<1, 5e-3, 5e-3)'
  [../]
[]

[Materials]

  ##materials for thermal evolution

  [q_plastic]
    type = ADComputePlasticWorkHeating
    beta_p = 0.5
    dirac_switch_react = dirac_switch_react
    thr_activation = 1
    use_PK2 = true
    use_lump = true
  []
  [q_elastic]
    type = ADComputeElasticWorkHeating
    temperature = temperature
    beta_av = 0.5
    dirac_switch_react = dirac_switch_react
    thr_activation = 1
    use_PK2 = true
    use_lump = true
  []

  #materials for switches
  [rate_dirac_switch_shock]
    type = ParsedMaterial
    f_name = rate_dirac_shock
    args = 'v_flag vx'
    function = 'if(v_flag > 0, 1., 0.)' 
  []
  [rate_dirac_switch_react]
    type = ParsedMaterial
    f_name = rate_dirac_react
    args = 'dirac_switch_shock v_flag time_react'
    function = 'if(dirac_switch_shock > 1, 1. / time_react, 0)'
  []
  [ADplate_const]
  	type = ADGenericConstantMaterial
  	prop_names = 'density specific_heat thermal_conductivity k_corrected'
  	prop_values = '1807e-3 2320e-6 0.67e-6 0.000159819' #kg/m3 kcal/mol-K W/m-K#
  []
  [elastic_tensor_plate]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []
  [compute_strain_plate]
    type = ComputeLagrangianStrain
    displacements = 'disp_x disp_y'
  []
  
  [./ComputeYdots]
    type = ADComputeYdotsRDXNS
    temperature = temperature
    Y1 = Y1
    Y2 = Y2

    Z1 = 15942.35
    Z2 = 2912.40

    E1 = 0.433477477
    E2 = 0.393899099

    Rg = 3.74505E-05
    T_trans = 1736.16
    a1 = -0.03203964
    b1 = 0.000591791

    a2 = 7.840288288
    b2 = -0.000722965
    
    dirac_switch_react = dirac_switch_react
    switch_react = 1 #when to start kinetics
    rate_limit = 1e10
    use_lump = true
    dynamic_tau = true
    thr_activation_rates = 1
  [../]
  
  [flow_stress_plate]
  	type = ADComputeLagrangianJCYieldStress
  	epsilon_ref = 1e3
  	temperature = temperature
  	T0 = 300
  	k = 3
  	compute = false
  	A = 0.1
  	B = 0.3
  	use_temp = 1.0
  	a_melt = 1.53
  	Tm0 = 552
  	n_h = 0.1
    melting_model = 0
    use_rate = 1.0
  []
  [compute_stress_plate]
    type = ADArtVisJ2Stress
    flow_stress_material = flow_stress_plate
    C0 = 0.1
    bulk = 1 #deprecated
    C1 = 1.0
    element_size = 1.9
    Yfinal = Y3
  []
  [compute_misternet_heat]
    type = ADComputeMISTERnetHeat
    dirac_switch_shock = dirac_switch_shock
    dirac_switch_react = dirac_switch_react
    temperature= temperature
    density = density
    specific_heat = specific_heat
    factor=1.0
    T_ref=300
    heat_time_react = 5 #ns
    heat_time_shock = 1 #ns
    vx = vx
    ax = ax
    thr_v = 0.5 #deprecated
    thr_a = 1 #deprecated
    direct_T = false
    dynamic_tau = true
    temp_crit = true
    use_sin = false
    element_size = 1.9
  []
  [pressureval]
    type = ADComputeIntPValRDXMISTERnetNSFull
    temperature = temperature
    element_size = 1.9
    T_ref = 300
    P0 = -6.02
    C0 = 0.1
    C1 = 1.0

    ##unreacted parameters
    A_u = 3174.2081
    B_u = -717.2443
    R1_u = 254.7515
    R2_u = 25.47515
    omega_u = 4.2306
    
    ##reacted parameters
    A_r = 923.8521
    B_r = -2375.7735
    R1_r = 11124.2648
    R2_r = 1112.42648
    omega_r = 2.3543

    #interpolation parameters
    Y_final = Y3
    use_fitted_eos = true
    flag_threshold = 1e-3
    
    use_RDX = 1.0
    
    ##HMX parameters
    Gamma = 0.7
    slope = 2.29
    A1 = 1668.9 #GPa
    A2 = 59.69 #GPa
    R1 = 5.9
    R2 = 2.1
    omega = 0.45
    
    ##call parameters
    thr_a = 0.1 ##acceleration smaller than this value CURRENTLY WORKING VALUE OF 2.5e-2
    thr_v = 0.5 ##velocity grater than this value
    vx = vx
    ax = ax

    #extrapolation flag
    extrapolate = 1.
    
    ## MISTERnet parameters
    up = 1.75 # hard coded for now
    density_i=density_i
    mask_size = 40
    use_mask = 0
	
	  ##CSV data
	  csv_shock = 'shock_new.csv'
	  csv_react = 'react_new.csv'
    csv_times = 'time_new.csv'
  []
[]

[AuxVariables]
  [./density_i]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [temperature_mister_shock]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = MaterialRealAux
      property = temperature_mister_shock
      variable = temperature_mister_shock
    []
  []
  [temperature_mister_react]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = MaterialRealAux
      property = temperature_mister_react
      variable = temperature_mister_react
    []
  []
  [time_react]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = MaterialRealAux
      property = time_react
      variable = time_react
    []
  []
  [heatrate_mister_shock]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADMaterialRealAux
      property = heatrate_mister_shock
      variable = heatrate_mister_shock
    []
  []

  [heatrate_mister_react]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADMaterialRealAux
      property = heatrate_mister_react
      variable = heatrate_mister_react
    []
  []
  [sxx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = cauchy_stress
      index_i = 0
      index_j = 0
    []
  []
  [sxy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = cauchy_stress
      index_i = 0
      index_j = 1
    []
  []
  [sxz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = cauchy_stress
      index_i = 0
      index_j = 2
    []
  []
  [syy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = cauchy_stress
      index_i = 1
      index_j = 1
    []
  []
  [syz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = cauchy_stress
      index_i = 1
      index_j = 2
    []
  []
  [szz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = cauchy_stress
      index_i = 2
      index_j = 2
    []
  []
  
  ##TEST FLOW DIRECTION PLOTTING
  
  [nxx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = flow_direction
      index_i = 0
      index_j = 0
    []
  []
  [nxy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = flow_direction
      index_i = 0
      index_j = 1
    []
  []
  [nxz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = flow_direction
      index_i = 0
      index_j = 2
    []
  []
  [nyy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = flow_direction
      index_i = 1
      index_j = 1
    []
  []
  [nyz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor =flow_direction
      index_i = 1
      index_j = 2
    []
  []
  [nzz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = flow_direction
      index_i = 2
      index_j = 2
    []
  []
  
  ##TEST CAUCHY FROM BLATZ-KO STRAIN ENERGY DENSITY FUNCTIONAL
  
  
  [exx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 0
      index_j = 0
    []
  []
  [exy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 0
      index_j = 1
    []
  []
  [exz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 0
      index_j = 2
    []
  []
  [eyy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 1
      index_j = 1
    []
  []
  [eyz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 1
      index_j = 2
    []
  []
  [ezz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 2
      index_j = 2
    []
  []
  

  ##elastic deformation gradient

  [Fexx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fe
      index_i = 0
      index_j = 0
    []
  []
  [Fexy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fe
      index_i = 0
      index_j = 1
    []
  []
  [Fexz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fe
      index_i = 0
      index_j = 2
    []
  []
  [Feyy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fe
      index_i = 1
      index_j = 1
    []
  []
  [Feyz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fe
      index_i = 1
      index_j = 2
    []
  []
  [Fezz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fe
      index_i = 2
      index_j = 2
    []
  []
  
  ##plastic deformation gradient

  [Fpxx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fp
      index_i = 0
      index_j = 0
    []
  []
  [Fpxy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fp
      index_i = 0
      index_j = 1
    []
  []
  [Fpxz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fp
      index_i = 0
      index_j = 2
    []
  []
  [Fpyy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fp
      index_i = 1
      index_j = 1
    []
  []
  [Fpyz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fp
      index_i = 1
      index_j = 2
    []
  []
  [Fpzz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Fp
      index_i = 2
      index_j = 2
    []
  []

  ##

  [Epxx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep
      index_i = 0
      index_j = 0
    []
  []
  [Epxy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep
      index_i = 0
      index_j = 1
    []
  []
  [Epxz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep
      index_i = 0
      index_j = 2
    []
  []
  [Epyy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep
      index_i = 1
      index_j = 1
    []
  []
  [Epyz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep
      index_i = 1
      index_j = 2
    []
  []
  [Epzz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep
      index_i = 2
      index_j = 2
    []
  []

  [Ep_dotxx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep_dot
      index_i = 0
      index_j = 0
    []
  []
  [Ep_dotxy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep_dot
      index_i = 0
      index_j = 1
    []
  []
  [Ep_dotxz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep_dot
      index_i = 0
      index_j = 2
    []
  []
  [Ep_dotyy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep_dot
      index_i = 1
      index_j = 1
    []
  []
  [Ep_dotyz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep_dot
      index_i = 1
      index_j = 2
    []
  []
  [Ep_dotzz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ep_dot
      index_i = 2
      index_j = 2
    []
  []

  [Eexx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee
      index_i = 0
      index_j = 0
    []
  []
  [Eexy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee
      index_i = 0
      index_j = 1
    []
  []
  [Eexz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee
      index_i = 0
      index_j = 2
    []
  []
  [Eeyy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee
      index_i = 1
      index_j = 1
    []
  []
  [Eeyz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee
      index_i = 1
      index_j = 2
    []
  []
  [Eezz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee
      index_i = 2
      index_j = 2
    []
  []

  [Ee_dotxx]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee_dot
      index_i = 0
      index_j = 0
    []
  []
  [Ee_dotxy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee_dot
      index_i = 0
      index_j = 1
    []
  []
  [Ee_dotxz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee_dot
      index_i = 0
      index_j = 2
    []
  []
  [Ee_dotyy]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee_dot
      index_i = 1
      index_j = 1
    []
  []
  [Ee_dotyz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee_dot
      index_i = 1
      index_j = 2
    []
  []
  [Ee_dotzz]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = RankTwoAux
      rank_two_tensor = Ee_dot
      index_i = 2
      index_j = 2
    []
  []

  ###jacobians

  [Je]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = RankTwoScalarAux
  		variable = Je
  		scalar_type = ThirdInvariant
      rank_two_tensor = Fe
      execute_on = timestep_end
  	[]
  []

  [detEe_dot]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = RankTwoScalarAux
  		variable = detEe_dot
  		scalar_type = FirstInvariant
      rank_two_tensor = Ee_dot
      execute_on = timestep_end
  	[]
  []

  [Jp]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = RankTwoScalarAux
  		variable = Jp
  		scalar_type = ThirdInvariant
      rank_two_tensor = Fp
      execute_on = timestep_end
  	[]
  []

  ##

  [pressure_av]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = pressure_av
      property = pressure_av
    []
  []
  
  ##TEST: JC variables output
  
  [theta]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = MaterialRealAux
      variable = theta
      property = theta
    []
  []
  [ep_rate]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = MaterialRealAux
      variable = ep_rate
      property = ep_rate
    []
  []
  [ep]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = MaterialRealAux
      variable = ep
      property = effective_plastic_strain
    []
  []
  ##output pressures
  [pressure_mg]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = pressure_mg
      property = pressure_mg
    []
  []
  [pressure_JWL]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = pressure_JWL
      property = pressure_JWL
    []
  []
  [pressure_total]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = pressure_total
      property = pressure_total
    []
  []
  [v_flag]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = MaterialRealAux
      variable = v_flag
      property = v_flag
    []
  []
  [r1]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = r1
      property = r1
    []
  []
  [r2]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = r2
      property = r2
    []
  []
  [q_decomposition]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = q_decomposition
      property = q_decomposition
    []
  []
  [q_elastic]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = q_elastic
      property = q_elastic
    []
  []
  [q_plastic]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = q_plastic
      property = q_plastic
    []
  []
  [Y1_dot_surrogate]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = Y1_dot_surrogate
      property = Y1_dot_surrogate
    []
  []
  [Y2_dot_surrogate]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = Y2_dot_surrogate
      property = Y2_dot_surrogate
    []
  []
  [Y3_dot_surrogate]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = Y3_dot_surrogate
      property = Y3_dot_surrogate
    []
  []
  [indicator_surrogate]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = MaterialRealAux
      variable = indicator_surrogate
      property = indicator_surrogate
    []
  []

  
  [us]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = ADMaterialRealAux
      variable = us
      property = us
    []
  []
  [called_up]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = MaterialRealAux
      variable = called_up
      property = called_up
    []
  []
  ##reaction Heats
  [Q1]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = ADMaterialRealAux
  		variable = Q1
  		property = Q1
  	[]
  []
  [Q2]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = ADMaterialRealAux
  		variable = Q2
  		property = Q2
  	[]
  []

  ##hotspot heat

  

  [ratedep]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = MaterialRealAux
  		variable = ratedep
  		property = ratedep
  	[]
  []
  [flow_stress]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = MaterialRealAux
  		variable = flow_stress
  		property = flow_stress
  	[]
  []

  ##

  [von_mises]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = RankTwoScalarAux
  		variable = von_mises
  		scalar_type = VonMisesStress
      rank_two_tensor = cauchy_stress
      execute_on = timestep_end
  	[]
  []
  [s_pressure]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = RankTwoScalarAux
  		variable = s_pressure
  		scalar_type = Hydrostatic
      rank_two_tensor = cauchy_stress
      execute_on = timestep_end
  	[]
  []
  [Y1dot_aux]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = TimeDerivativeAux
      variable = Y1dot_aux
      functor = Y1
    []
  []
  [Y2dot_aux]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = TimeDerivativeAux
      variable = Y2dot_aux
      functor = Y2
    []
  []
  [Y3dot_aux]
    family = MONOMIAL
    order = CONSTANT
    [AuxKernel]
      type = TimeDerivativeAux
      variable = Y3dot_aux
      functor = Y3
    []
  []
[]

[Executioner]
	type = Transient
	#dt = 5e-3
  [./TimeStepper]
    type = FunctionDT
    function = timestep
    min_dt = 1e-8
  [../]
	nl_rel_tol = 1E-6
	nl_abs_tol = 1E-6
	solve_type = Newton
	#timestep_tolerance = 1e-6
	#[TimeIntegrator]
	#	type = NewmarkBeta
	#	inactive_tsteps = 5
	#[]
	petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter -snes_type'
	petsc_options_value = '401 hypre boomeramg 60 vinewtonrsls' 
	automatic_scaling = true
	line_search = 'none'
[]

#[Preconditioning]
#	[full]
#		type = SMP
#		full = true
#	[]
#[]


[Outputs]
  exodus = true
  interval = 20
[]