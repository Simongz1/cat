E = 0.87e-3
nu = 0.45
width = 4

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = true
  block = 1
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [temperature]
  []
  [c]
    #initial_condition = 0
  []
[]

[AuxVariables]
	##dynamics variables

	[./vx]
	[../]
	[./ax]
	[../]
	[./vy]
	[../]
	[./ay]
	[../]
  [./vz]
	[../]
	[./az]
	[../]

  [./gc]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dummyc]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Bounds]
  [c_high]
    type = ConstantBounds
    variable = dummyc
    bounded_variable = c
    bound_type = upper
    bound_value = 1
  []
  [c_old]
    type = VariableOldValueBounds
    variable = dummyc
    bounded_variable = c
    bound_type = lower
  []
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = 'circle.msh'
  []
[]

#[Mesh]
#    type = GeneratedMesh
#    dim = 3
#    nx = 35
#    ny = 200
#    nz = 1
#    xmax = 3
#    xmin = 0
#    ymax = 20
#    ymin = 0
#    zmax = 0.2
#    zmin = 0
#    boundary_id = '0 1 2 3 4 5'
#    boundary_name = 'back bottom right top left front'
#    add_subdomain_ids = 2
#[]

[MeshModifiers]
  [moving_circle]
    type = CoupledVarThresholdElementSubdomainModifier
    coupled_var = c
    criterion_type = ABOVE
    threshold = 0.95
    subdomain_id = 2
    execute_on = timestep_end
    block = 1
    reinitialize_subdomains = 1 #this only reinitializes the alive block
  []
[]

[Kernels]
  ##temperature kernels
  [dTdt]
  	type = ADHeatConductionTimeDerivative
    variable = temperature
    density_name = 'density'
    specific_heat = 'specific_heat'
  []
  [nabla2T]
  	type = ADHeatConduction
    variable = temperature
  []
  [FullHS]
    type = LIPITHS
    variable = temperature
    beta_p = 0.5
    beta_av = 0.5
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
  

  #balance of linear momentum
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


  ##fracture kernels
  [c_dot]
    type = ADTimeDerivative
    variable = c
  []
  [gc_l_gradc]
    type = ADMatDiffusion
    variable = c
    diffusivity = gcl
  []
  [gcc]
    type = ADMatBodyForce
    variable = c
    material_property = gcc
  []
  [dDdcH]
    type = ADMatBodyForce
    variable = c
    material_property = dDdcH
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
  [./velocity]
  	type = PresetVelocity
  	boundary = left
  	velocity = -1
  	variable = disp_x
  [../]
  #[plain]
  #  type = DirichletBC
  #  variable = disp_z
  #  boundary = 'front back'
  #  value = 0
  #[]
  [plain2]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0
  []
  #[plane_right]
  #  type = DirichletBC
  #  variable = disp_y
  #  boundary = 'front back'
  #  value = 0
  #[]
  [fix]
    type = DirichletBC
    variable = disp_y
    boundary = 'left right'
    value = 0
  []
[]

[ICs]
	[temp_IC]
		type = ConstantIC
		variable = temperature
		value = 300
	[]
  [ConstantIC]
    type = ConstantIC
    variable = gc
    value = 100e-6
  []
  #[thumb]
  #  type = ThumbIC
  #  height = 5
  #  invalue = 1
  #  outvalue = 0
  #  width = 0.5
  #  xcoord = 1.5
  #  variable = c
  #[]
[]

[Materials]
  ##materials on block 1: plate
  [plate_const]
  	type = ADGenericConstantMaterial
  	prop_names = 'density specific_heat thermal_conductivity alpha'
  	prop_values = '1100e-9 1250 0.13e-6 0.5e-4' #kg/m3 J/kg-K W/m-K#
  []
  [other_constants]
  	type = ADGenericConstantMaterial
  	prop_names = 'visco kdamage nu'
  	prop_values = '0.1 1e-6 ${nu}'
  []
  [l]
    type = ADParsedMaterial
    property_name = l
    coupled_variables = 'h_min'
    expression = '${width} * h_min'
  []
  [elastic_tensor_plate]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []
  [compute_strain_plate]
    type = ADComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
    #outputs = exodus
  []
  [flow_stress]
    type = ADDerivativeParsedMaterial
    property_name = flow_stress
    constant_names = 'A B n m'
    constant_expressions = '0.3e-4 0.5e-4 1 1'
    material_property_names = 'ep theta(temperature)'
    coupled_variables = 'c temperature'
    expression = '(A + B * ep ^ n) * (1 - theta ^ m)'
    compute = false
    #outputs = exodus
  []
  [theta]
    type = ADDerivativeParsedMaterial
    property_name = theta
    constant_names = 'T0 Tmelt'
    constant_expressions = '300 500'
    coupled_variables = 'temperature'
    expression = 'min(0.9, (temperature - T0)/(Tmelt - T0))'
    #outputs = exodus
  []
  [compute_stress_plate]
    type = ADComputeChainEnergy
    flow_stress_material = flow_stress
    C0 = 0.1
    C1 = 1.0
    c = c
    gc = gc
    displacements = 'disp_x disp_y disp_z'
    h_min = h_min
    outputs = exodus
  []
  
  ##elastic energy material
  [degradation_material]
    type = ADDerivativeParsedMaterial
    material_property_names = 'kdamage'
    coupled_variables = 'c'
    expression = '(1 - c) ^ 2 + kdamage'
    derivative_order = 1
    property_name = D
    outputs = exodus
  []
  [dDdcH]
    type = ADDerivativeParsedMaterial
    property_name = dDdcH
    coupled_variables = 'c'
    material_property_names = 'dDdc:=D[D,c] visco Hist'
    expression = '- dDdc * Hist  / visco'
    derivative_order = 1
    outputs = exodus
  []
  [gcc]
    type = ADDerivativeParsedMaterial
    property_name = gcc
    coupled_variables = 'c gc'
    material_property_names = 'visco l'
    expression = '- gc * c / (l * visco)'
    derivative_order = 1
  []
  [gcl]
    type = ADDerivativeParsedMaterial
    property_name = gcl
    material_property_names = 'l visco'
    coupled_variables = 'gc'
    expression = 'gc * l / visco'
    derivative_order = 1
  []

  ##test
  [dDdc_out]
    type = ADDerivativeParsedMaterial
    property_name = dDdc_test
    coupled_variables = 'c'
    material_property_names = 'dDdc:=D[D,c]'
    expression = 'dDdc'
    outputs = exodus
  []
[]

[AuxKernels]
  [./min]
    type = ElementLengthAux
    variable = h_min
    method = min
    execute_on = 'initial timestep_end'
  [../]
  [./max]
    type = ElementLengthAux
    variable = h_max
    method = max
    execute_on = 'initial timestep_end'
  [../]
  
[]

[AuxVariables]
  [h_min]
    order = CONSTANT
    family = MONOMIAL
  []
  [h_max]
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
[]

[Executioner]
  type = Transient
  line_search = none
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter -snes_type'
  petsc_options_value = '401 hypre boomeramg 60 vinewtonrsls' 
  automatic_scaling = true
  solve_type = Newton
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-7
  start_time = 0.0
  dt = 4e-3
  end_time = 1e10
[]

[Postprocessors]
  [./react_left]
    type = ADSidesetReaction
    direction = '-1 0 0'
    stress_tensor = stress
    boundary = left
  [../]
  [./react_right]
    type = ADSidesetReaction
    direction = '1 0 0'
    stress_tensor = stress
    boundary = right
  [../]
[]

[Outputs]
  exodus = true
  time_step_interval = 10
[]


