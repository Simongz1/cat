#try to keep these fixed
E = 20.76
nu = 0.24

{UNREACTED_DISTRIBUTIONS}

{REACTED_DISTRIBUTIONS}

[GlobalParams]
  #general
  displacements = 'disp_x disp_y disp_z'

  #this should never change
  bulk_MicroID = 99

  #for microstructure
  bulk_grains = {BULK_GRAINS}
  euler_angles = false

  dirac_tolerance = 0
  use_mixture = {USE_MIXTURE}
  
  #define binder mechanical properties
  binder_bulk = {BINDER_BULK}
  binder_yield = {BINDER_YIELD}
  binder_shear = {BINDER_SHEAR}

  fraction_csv = fraction_csv
  use_tabular_time = {TABULAR_TIME}
  use_distributions = {DISTRIBUTIONS}
  consistent_tau = true

  use_loaded_microstructure = {LOADED_MICROSTRUCTURE}
  use_complete_burn = {COMPLETE_BURN}

  #testing bulk grains
  n_cracks = 1 #unused
  l_cracks = 1 #unused
  n_pores = 10000000
  pore_probability = {POROSITY} # default is 0.01

  bulk_RDX_fraction = 1. ##saturated RDX
  pore_RDX_fraction = 1. ##this should be deprecated
  range_pore = '101 103' #keep constant

  ##csv data assignment
  csv_shock_pore = 'shock_pore.csv'
  csv_react_pore = 'react_pore.csv'
  csv_times_pore = 'time_pore.csv'
  csv_density_pore = 'density_pore.csv'
  csv_fraction_pore = 'fraction_pore.csv'

  #use old time formulation
  use_gating = {USE_GATING} 

  #scale reaction time
  tau_react_scaling = {SCALING_TAU}
[]

#add a block for the loaded microstructure in case it is needed

{MICROSTRUCTURE_FUNCTION}

{MICROSTRUCTURE_VARIABLE}

{MICROSTRUCTURE_IC}

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
  [disp_z]
    order = FIRST
  []

  [temperature]
    order = FIRST
  []
  [tracking]
    order = CONSTANT
    family = MONOMIAL
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
  [./vz]
    order = FIRST
    family = LAGRANGE
	[../]
	[./az]
    order = FIRST
    family = LAGRANGE
	[../]
[]

[Mesh]
  type = GeneratedMesh
  dim = 3 #this is fixed
  nx = {ELEM_PERP_1}
  ny = {ELEM_SHOCK_DIR}
  nz = {ELEM_PERP_2}
  xmax = {DIM_PERP_1}
  xmin = 0
  ymax = {DIM_SHOCK_DIR}
  ymin = 0
  zmax = {DIM_PERP_2}
  zmin = {ZMIN}
  boundary_id = '0 1 2 3 4 5' #keep these values
  boundary_name = 'back bottom right top left front' #keep these values
  elem_type = HEX8
[]

[Kernels]
  [tracking_t]
    type = ADTimeDerivative
    variable = tracking
  []

  [dirac_switch_shock_time_derivative]
    type = ADTimeDerivative
    variable = dirac_switch_shock
  []
  [dirac_switch_react_time_derivative]
    type = ADTimeDerivative
    variable = dirac_switch_react
  []
  [rate_shock]
    type = ADMatBodyForce
    variable = dirac_switch_shock
    material_property = rate_dirac_shock
  []
  [rate_react]
    type = ADMatBodyForce
    variable = dirac_switch_react
    material_property = rate_dirac_react
  []
  [rate_tracking]
    type = ADMatBodyForce
    variable = tracking
    material_property = rate_tracking
  []
	
  ##temperature kernels
  [dTdt]
  	type = MassLumpedTimeDerivative
    variable = temperature
    #specific_heat = 'specific_heat'
  []
  [nabla2T]
  	type = ADHeatConduction
    variable = temperature
    thermal_conductivity = 'k_corrected'
  []

  ##MISTERNET HEATS
  [heat_shock]
    type = ADMatBodyForce
    variable = temperature
    material_property = scaled_shock
  []
  [heat_react]
    type = ADMatBodyForce
    variable = temperature
    material_property = scaled_react
  []
  [sdx]
    type = ADDynamicStressDivergenceTensors
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
  []
  [sdy]
    type = ADDynamicStressDivergenceTensors
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
  []
  [sdz]
    type = ADDynamicStressDivergenceTensors
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
  []

	[./inertia_x]
		type = ADInertialForce
		variable = disp_x
		velocity = vx
		acceleration = ax
		beta = 0.3025 
		gamma = 0.6 
	[../]
	[./inertia_y]
		type = ADInertialForce
		variable = disp_y
		velocity = vy
		acceleration = ay
		beta = 0.3025 
		gamma = 0.6 
	[../]
  [./inertia_z]
		type = ADInertialForce
		variable = disp_z
		velocity = vz
		acceleration = az
		beta = 0.3025 
		gamma = 0.6 
	[../]

  [chem_from_prop]
    type = ADMatHeatSource
    variable = temperature
    material_property = scaled_chem
    scalar = 1 #can switch sign
  []
  [PlasticHeat]
    type = ADMatHeatSource
    variable = temperature
    material_property = scaled_plastic
    scalar = 1 #to assign sign
  []
  [PressureHeat]
    type = ADMatHeatSource
    variable = temperature
    material_property = scaled_elastic
    scalar = 1
  []
	[./dY1_dt]
		type = MassLumpedTimeDerivative
		variable = Y1
	[../]
	[./dY2_dt]
		type = MassLumpedTimeDerivative
		variable = Y2
	[../]
	[./dY3_dt]
		type = MassLumpedTimeDerivative
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
  [Y1_surr]
    type = ADMatBodyForce
    variable = Y1
    material_property = Y1_dot_surrogate
  []
  [Y2_surr]
    type = ADMatBodyForce
    variable = Y2
    material_property = Y2_dot_surrogate
  []
  [Y3_surr]
    type = ADMatBodyForce
    variable = Y3
    material_property = Y3_dot_surrogate
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
  [./vz]
    type = NewmarkVelAux
    variable = vz
    acceleration = az
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
  [./az]
    type = NewmarkAccelAux
    variable = az
    displacement = disp_z
    velocity = vz
    beta = 0.3025
	[../]	
[]

[BCs]
  [x_fix]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  []
  [z_fix]
    type = DirichletBC
    variable = disp_z
    boundary = 'front back'
    value = 0.0
  []
  [./velocity]
    type = PresetVelocity
    boundary = bottom
    velocity = {IMPACT_VELOCITY}
    variable = disp_y
  [../]
  [./top]
    type = DirichletBC
    boundary = top
    value = 0
    variable = disp_y
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
[]

[Materials]

  [scaled_shock]
    type = ADParsedMaterial
    material_property_names = 'heatrate_mister_shock density specific_heat'
    expression = 'heatrate_mister_shock / (density * specific_heat)'
    property_name = scaled_shock
  []
  [scaled_react]
    type = ADParsedMaterial
    material_property_names = 'heatrate_mister_react density specific_heat'
    expression = 'heatrate_mister_react / (density * specific_heat)'
    property_name = scaled_react
  []
  [scaled_elastic]
    type = ADParsedMaterial
    material_property_names = 'q_elastic density specific_heat'
    expression = 'q_elastic / (density * specific_heat)'
    property_name = scaled_elastic
  []
  [scaled_plastic]
    type = ADParsedMaterial
    material_property_names = 'q_plastic density specific_heat'
    expression = 'q_plastic / (density * specific_heat)'
    property_name = scaled_plastic
  []
  [scaled_chem]
    type = ADParsedMaterial
    material_property_names = 'q_decomposition density specific_heat'
    expression = 'q_decomposition/ (density * specific_heat)'
    property_name = scaled_chem
  []
  [artificial_diff_temp]
    type = ADParsedMaterial
    material_property_names = 'norm_gradT thermal_conductivity'
    expression = 'thermal_conductivity * (1 + 5e-1*norm_gradT)'
    property_name = 'D'
    
  []
  ########DEFINE EQUATIONS OF STATE HERE
  [dPdT]
    type = ADDerivativeParsedMaterial
    material_property_names = 'density specific_heat Je time_react'
    coupled_variables = 'temperature Y1 dirac_switch_react'
    expression = 'if(dirac_switch_react > time_react, (density * specific_heat / Je) * (Y1 * omega_unreacted + (1 - Y1) * omega_reacted), 0.0)'
    derivative_order = 1
    constant_names = 'omega_unreacted omega_reacted'
    constant_expressions = '0.37 0.77' #keep constant for now
    property_name = dPdT
    
  []
  
  [q_plastic]
    type = ADComputePlasticWorkHeating
    beta_p = 0.5
    dirac_switch_react = dirac_switch_react
    thr_activation = 1
    
  []
  [q_elastic]
    type = ADComputeElasticWorkHeating
    temperature = temperature
    beta_av = 0.5
    dirac_switch_react = dirac_switch_react
    thr_activation = 1
    
  []

  #materials for switches
  [rate_dirac_switch_shock]
    type = ADDerivativeParsedMaterial
    property_name = rate_dirac_shock
    material_property_names = 'v_flag time_shock'
    function = 'if(v_flag > 0, 1., 0.)'
  []
  [rate_dirac_switch_react]
    type = ADDerivativeParsedMaterial
    property_name = rate_dirac_react
    coupled_variables = 'dirac_switch_shock'
    material_property_names = 'time_react time_shock'
    function = 'if(dirac_switch_shock >= time_shock, 1., 0.)'
  []
  [ADplate_const]
  	type = ADGenericConstantMaterial
  	prop_names = 'specific_heat thermal_conductivity'
  	prop_values = '2320e-6 0.37e-6' 
  []
  [k_corrected]
    type = ADDerivativeParsedMaterial
    property_name = k_corrected
    material_property_names = 'specific_heat thermal_conductivity density'
    expression = 'thermal_conductivity / (density * specific_heat)' 
  []
  [elastic_tensor_plate]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []
  [compute_strain_plate]
    type = ADComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
  []

  [yield_mixture]
    type = ADDerivativeParsedMaterial
    derivative_order = 2
    expression = 'fraction_csv * (A + B * ep ^ (1+n)) * (1 - theta^m) + (1 - fraction_csv) * yield_binder'
    material_property_names = 'ep theta(temperature)'
    constant_names = 'A B n m yield_binder'
    constant_expressions = '0.3 0.1 0.1 3 {BINDER_YIELD}'
    property_name = yield_mixture
    coupled_variables = 'temperature fraction_csv'
    compute = false
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
    thr_activation_rates = 0
  [../]

  [compute_stress_plate]
    type = ADPBXStress
    flow_stress_material = yield_mixture
    poisson_binder = {BINDER_POISSON}
    C0 = 0.1
    C1 = 1.0
    element_size = {SHOCKDIR}
    Yinitial = Y1
    displacements = 'disp_x disp_y disp_z'

    A_unreacted = 13398.42268091
    R1_unreacted = 68.1217746
    R2_unreacted = 6.81217746
    B_unreacted = 3810.3712782
    omega_unreacted = 0.37

    #reacted
    A_reacted = 1000
    R1_reacted = 59.99058313
    R2_reacted = 5.999058313
    B_reacted = 2738.01645291
    omega_reacted = 0.77

    temperature = temperature
    {OUTPUTS}
  []
  [compute_misternet_heat]
    type = ADComputeMISTERnetHeat
    dirac_switch_shock = dirac_switch_shock
    dirac_switch_react = dirac_switch_react
    temperature = temperature
    density = density
    specific_heat = specific_heat
    T_ref = 300
    element_size = {SHOCKDIR}
    fraction_csv = fraction_csv

    ##for velocity and acceleration thresholds
    v_components = 'vx vy vz'
    a_components = 'ax ay az'
    {OUTPUTS}
  []
  [pressureval]
    type = ADComputeIntPValRDXMISTERnetNSFull
    element_size = {SHOCKDIR}
    tracking = tracking    
    ##call parameters
    thr_a = 1.5 ##acceleration smaller than this value CURRENTLY WORKING VALUE OF 2.5e-2
    thr_v = 0.5 ##velocity grater than this value
    v_components = 'vx vy vz'
    a_components = 'ax ay az'

    density_i = density_i	

	  ###########CSV data
	  csv_shock = 'shock_new.csv'
	  csv_react = 'react_new.csv'
    csv_times = 'time_new.csv'

    ##variable that stores density from csv
    density_csv = density_csv

    #############
    csv_unreacted = 'unreacted.csv'
    csv_reacted = 'reacted.csv'
    {OUTPUTS}
  []
[]

[AuxKernels]
  #######to assign density as a function of the voronoi tessellation
  [readVoronoiStructure]
    type = PolycrystalDensityAux
    variable = density_csv
    execute_on = INITIAL
    csv_density = 'densities_new.csv'
    density_i = density_i
    density_scaling = 1e-3

    ##testing bulk grains
  
    bulk_RDX_density = {BULK_DENSITY} #provided in g/cm³
  []
  [readVoronoiRDXFraction]
    type = PolycrystalFractionAux
    variable = fraction_csv
    execute_on = 'INITIAL TIMESTEP_END'
    csv_fraction = 'fractions_new.csv'
    density_i = density_i  
    bulk_RDX_fraction = 1.
  []
[]

[UserObjects]
  [voronoi_density]
    type = PolycrystalDensityUO
    num_grains = 800 #keep this number high enough
    target_grains = '1'
    range_in = '{PARTICLE_RANGE_SMALL} {PARTICLE_RANGE_LARGE}' #this is for particles
    range_out = '{BINDER_RANGE_SMALL} {BINDER_RANGE_LARGE}' #this is for binder
    generate_matrix = true #keep true
    max_grain_size = 1000
    execute_on = 'INITIAL TIMESTEP_END'
    min_center_spacing = {SMALL_SIZE}
    matrix_thickness = {BINDER_WIDTH}
    sizes = '{SMALL_SIZE} {BIG_SIZE}'
    sizes_fraction = '{SMALL_PROPORTION} {BIG_PROPORTION}'

    ##for fraction assignment
    csv_fraction = 'fractions_new.csv'
  []
[]

##set global params for the pore assignment
[GlobalParams]
  
[]

[AuxVariables]
  [./density_i]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [density_csv]
    order = CONSTANT
    family = MONOMIAL
  []
  [fraction_csv]
    order = CONSTANT
    family = MONOMIAL
  []
  [grainID]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler1]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler2]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler3]
    order = CONSTANT
    family = MONOMIAL
  []

  [von_mises]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = ADRankTwoScalarAux
  		variable = von_mises
  		scalar_type = VonMisesStress
      rank_two_tensor = stress
      execute_on = timestep_end
  	[]
  []
  [s_pressure]
  	family = MONOMIAL
  	order = CONSTANT
  	[AuxKernel]
  		type = ADRankTwoScalarAux
  		variable = s_pressure
  		scalar_type = Hydrostatic
      rank_two_tensor = stress
      execute_on = timestep_end
  	[]
  []
[]

[Executioner]
	type = Transient
  [./TimeStepper]
    type = FunctionDT
    function = timestep
    min_dt = 1e-8
  [../]
	nl_rel_tol = 1E-8
	nl_abs_tol = 1E-8
	solve_type = Newton
	petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter -snes_type'
	petsc_options_value = '401 hypre boomeramg 60 vinewtonrsls' 
	automatic_scaling = true
	line_search = 'none'
[]

[Functions]
  [./timestep]
    type = ParsedFunction
    value = 'if(t<1, 5e-3, 5e-3)'
  [../]
[]

[Outputs]
  exodus = true
  time_step_interval = 20
  [./mycheckpoints]
    type = Checkpoint
    wall_time_interval = {{FREC}}
  []
[]