[GlobalParams]

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
  [temperature]
    order = FIRST
  []

  [tau_hotspot]
    order = FIRST
    family = LAGRANGE
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

[]

[Mesh]
	type = GeneratedMesh
	dim = 2
	nx = 1
	ny = 1
	xmax = 1.9 #nm
	xmin = 0
	ymax = 1.1 #nm
	ymin = 0
	boundary_id = '0 1 2 3'
	boundary_name = 'bottom right top left'
  elem_type = QUAD4
[]

[Kernels]

  [./d_dirac_switch_dt_shock]
    type = ADTimeDerivative
    variable = dirac_switch_shock
  [../]
  [./d_dirac_switch_dt_react]
    type = ADTimeDerivative
    variable = dirac_switch_react
  [../]
  [./dirac_switch_shock]
		type = MaskedBodyForce
		variable = dirac_switch_shock
		value = 1.0 
		function = dirac_switch_induction_time_shock
		mask = dirac_switch_pressure_shock
	[../]
	[./dirac_switch_react]
		type = MaskedBodyForce
		variable = dirac_switch_react
		value = 1.0 
		function = dirac_switch_induction_time_react
		mask = dirac_switch_pressure_react
	[../]
	
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

  ##heat sources for plastic dissipation, elastic dissipation, and chemical decomposition
  [chem_from_prop]
    type = ADMatHeatSource
    variable = temperature
    material_property = q_decomposition
    scalar = 1 #can switch sign
  []

  [tau_hotspot]
    type = ADTimeDerivative
    variable = tau_hotspot
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

[]

[AuxKernels]
	
[]

[BCs]

[]

[ICs]
  [constantT]
    variable = temperature
    type = ConstantIC
    value = 1055.066986
  []
  [tau_hotspot]
    type = ConstantIC
    variable = tau_hotspot
    value = 1
  []

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
[]

[Functions]
  [./dirac_switch_induction_time_shock]
    type = ParsedFunction
    value = '1. / 1.' # 1/ns
  [../]
  [./dirac_switch_induction_time_react]
    type = ParsedFunction
    value = '1. / 5.' # 1/ns
  [../]
  [./timestep]
    type = ParsedFunction
    value = 'if(t<1, 1e-5, 1e-5)'
  [../]
[]

[Materials]
  [./w_time_dirac_switch_shock]
    type = ParsedMaterial
    f_name = dirac_switch_pressure_shock
    function = '1' 
  [../]
  [./w_time_dirac_switch_react]
    type = ParsedMaterial
    f_name = dirac_switch_pressure_react
    function = '1'
  [../]
  [ADplate_const]
  	type = ADGenericConstantMaterial
  	prop_names = 'density specific_heat thermal_conductivity k_corrected'
  	prop_values = '1807e-3 2320e-6 0.67e-6 0.000159819' #kg/m3 kcal/mol-K W/m-K#
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
    switch_react = 0 #when to start kinetics
    rate_limit = 1e10
    use_lump = true
  [../]
[]

[AuxVariables]
  [./density_i]
    order = CONSTANT
    family = MONOMIAL
  [../]

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
  [./TimeStepper]
    type = FunctionDT
    function = timestep
    min_dt = 1e-8
  [../]
	nl_rel_tol = 1E-6
	nl_abs_tol = 1E-6
	solve_type = Newton
	petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter -snes_type'
	petsc_options_value = '401 hypre boomeramg 60 vinewtonrsls' 
	automatic_scaling = true
	line_search = 'none'
[]

[Preconditioning]
	[full]
		type = SMP
		full = true
	[]
[]


[Outputs]
  exodus = true
  interval = 20
[]