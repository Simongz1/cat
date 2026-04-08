import pandas as pd
import numpy as np
import subprocess as su
from pathlib import Path

#provide dataframe with all the data for distributions
df = pd.read_csv('distributions.csv')

#define simulation dimension
dim = (input('Provide simulation dimension (2D) or (3D): '))
vel = float(input('Provide impact velocity (km/s): '))
elemsize = np.array((input('Voxel dimensions in microns for the shock and perpendicular directions: ')).split(','), dtype = float)

#query PBX or random
pbx = input('Type of microstructure to generate (PBX), (RANDOM) or (LOADED): ')

match pbx:
    case 'PBX':
        particle_sizes = np.array((input('Provide desired particle sizes (microns) for a bimodal distribution: ')).split(','), dtype = float)
        particle_proportions = np.array((input('Provide fractions of small and big particle sizes: ')).split(','), dtype = float)
        particle_poro = float(input('Provide particle porosity (%): ')) / 100
        binder_width = (input('Provide desired binder width (microns): '))
        binder_range = np.array((input('Provide minimum and maximum microstructure value for binder (0, 92): ')).split(','), dtype = int)
        particle_range = binder_range

        #define the rest of paramters
        bulk_grains = str('true')
        use_mixture = input('Use mixture for mechanics (true) or (false): ')
        binder_properties = np.array((input('Provide binder bulk, yield, shear modulus (GPa), and poisson modulus: ')).split(','), dtype = float)
        
        tabular_time = input('Use tabular time distribution? (true) or (false): ') #default to value
        use_distributions = input('Use gaussian time distributions: (true) or (false): ')
        use_gating = input('Use gating for heat and chemical sources (true) or (false): ')
        loaded_microstructure = str('false')
        datafile = 'a'

        #define block for function and variable
        block_function = """
            [Functions]
                [loaded_microstructure]
                    type = ParsedFunction
                    expression = a
                    symbol_names = 'a'
                    symbol_values = '1'
                []
            []

        """
        block_variable = """
            [AuxVariables]
                [loaded_microstructure]
                    order = CONSTANT
                    family = MONOMIAL
                []
            []
        """
        block_IC = """

        """
        
    case 'RANDOM':
        range_microstructures = np.array((input('Provide minimum and maximum microstructure ID: ').split(',')), dtype = int)
        particle_sizes = [50, 50] #default
        particle_proportions = [0.5, 0.5] #default
        particle_poro = 0.01 #default
        binder_width = 5 #default
        binder_range = range_microstructures
        particle_range = range_microstructures

        bulk_grains = str('false') #default

        use_mixture = input('Use mixture for mechanics (true) or (false): ')
        binder_properties = np.array((input('Provide binder bulk, yield, shear modulus (GPa), and poisson modulus: ')).split(','), dtype = float)
        tabular_time = input('Use tabular time distribution? (true) or (false): ') #default to value
        use_distributions = input('Use gaussian time distributions? (true) or (false): ')
        use_gating = input('Use gating for heat and chemical sources? (true) or (false): ')
        loaded_microstructure = str('false')
        datafile = 'a'
        
        #define block for function and variable
        block_function = """
            [Functions]
                [loaded_microstructure]
                    type = ParsedFunction
                    expression = a
                    symbol_names = 'a'
                    symbol_values = '1'
                []
            []
        """
        block_variable = """
            [AuxVariables]
                [loaded_microstructure]
                    order = CONSTANT
                    family = MONOMIAL
                []
            []
        """
        block_IC = """

        """
    
    case 'LOADED':
        datafile = input('DATA.txt file name without extension: ')
        if (dim == '2D'):
            zloc = float(input('Provide z location of the slice to extract: '))
        binder_range = np.array((input('Provide minimum and maximum microstructure value for binder (0, 92): ')).split(','), dtype = int)
        range_microstructures = binder_range
        bulk_grains = str('true')
        loaded_microstructure = str('true')

        particle_sizes = [50, 50] #default
        particle_proportions = [0.5, 0.5] #default
        particle_poro = 0.01 #default
        binder_width = 5 #default
        particle_range = binder_range

        use_mixture = input('Use mixture for mechanics (true) or (false): ')
        binder_properties = np.array((input('Provide binder bulk, yield, shear modulus (GPa), and poisson modulus: ')).split(','), dtype = float)

        tabular_time = input('Use tabular time distribution? (true) or (false): ') #default to value
        use_distributions = input('Use gaussian time distributions: (true) or (false): ')
        use_gating = input('Use gating for heat and chemical sources (true) or (false): ')

        #define block for function and variable
        block_function = """
            [Functions]
                [loaded_microstructure]
                    type = PiecewiseMultilinear
                    data_file = '{datafile}.txt'
                []  
            []
        """
        block_variable = """
            [AuxVariables]
                [loaded_microstructure]
                    order = CONSTANT
                    family = MONOMIAL
                []
            []
        """
        block_IC = """
            [ICs]
                [loaded_microstructure]
                    type = FunctionIC
                    function = loaded_microstructure
                    variable = loaded_microstructure
                []
            []
        """

#request complete burn
complete_burn = input('Use complete burn model? (true) or (false): ')

#template block for unreacted
block_temp_unreacted = """  
    [{name}_unreacted] type = Normal mean = {mean} standard_deviation = {std_dev} []
"""
#template block for reacted
block_temp_reacted = """  
    [{name}_reacted] type = Normal mean = {mean} standard_deviation = {std_dev} []
"""
blocks_unreacted = []
blocks_reacted = []

#iterate over values to form templated blocks
for _, row in df.iterrows():
    block_unreacted = block_temp_unreacted.format(name=int(row['Up']*10),
                              mean=row['mean_nonreactive'] if pd.notna(row['mean_nonreactive']) else 0,
                              std_dev=row['std_nonreactive'] if pd.notna(row['std_nonreactive']) else 1e-8)
    blocks_unreacted.append(block_unreacted)

for _, row in df.iterrows():
    block_reacted = block_temp_reacted.format(name=int(row['Up']*10),
                                              mean=row['mean_reactive'] if pd.notna(row['mean_reactive']) else 0,
                                              std_dev=row['std_reactive'] if pd.notna(row['std_reactive']) else 1e-8)
    blocks_reacted.append(block_reacted)

#format all the blocks
unreacted_out = "[Distributions]\n" + "".join(blocks_unreacted) + "[]\n"
reacted_out = "[Distributions]\n" + "".join(blocks_reacted) + "[]\n"

#write into file
with open("distributions_template.i", "r") as f:
    template_text = f.read()

#append new text
final = template_text.replace("{UNREACTED_DISTRIBUTIONS}", unreacted_out)
final = final.replace("{REACTED_DISTRIBUTIONS}", reacted_out)

#replace here what is common for both cases
final = final.replace("{IMPACT_VELOCITY}", str(vel))

#decompose binder properties
binder_bulk, binder_yield, binder_shear, binder_poisson = binder_properties

final = final.replace("{BINDER_BULK}", str(binder_bulk))
final = final.replace("{BINDER_YIELD}", str(binder_yield))
final = final.replace("{BINDER_SHEAR}", str(binder_shear))
final = final.replace("{BINDER_POISSON}", str(binder_poisson))

#replace the rest here
final = final.replace("{BULK_GRAINS}", bulk_grains)
final = final.replace("{USE_MIXTURE}", str(use_mixture))
final = final.replace("{TABULAR_TIME}", tabular_time)
final = final.replace("{DISTRIBUTIONS}", str(use_distributions))
final = final.replace("{PARTICLE_RANGE_SMALL}", str(particle_range[0]))
final = final.replace("{PARTICLE_RANGE_LARGE}", str(particle_range[-1]))
final = final.replace("{BINDER_RANGE_SMALL}", str(binder_range[0]))
final = final.replace("{BINDER_RANGE_LARGE}", str(binder_range[-1]))
final = final.replace("{SMALL_SIZE}", str(particle_sizes[0]))
final = final.replace("{BIG_SIZE}", str(particle_sizes[-1]))
final = final.replace("{BINDER_WIDTH}", str(binder_width))
final = final.replace("{SMALL_PROPORTION}", str(particle_proportions[0]))
final = final.replace("{BIG_PROPORTION}", str(particle_proportions[-1]))
final = final.replace("{POROSITY}", str(particle_poro))
final = final.replace("{USE_GATING}", str(use_gating))

#new
final = final.replace("{LOADED_MICROSTRUCTURE}", str(loaded_microstructure))
final = final.replace("{COMPLETE_BURN}", str(complete_burn))
final = final.replace("{MICROSTRUCTURE_FUNCTION}", block_function.format(datafile=str(datafile)))
final = final.replace("{MICROSTRUCTURE_VARIABLE}", block_variable)
final = final.replace("{MICROSTRUCTURE_IC}", block_IC)
final = final.replace("{SHOCKDIR}", str(elemsize[0]))

#for mesh generation
match dim:
    case '2D':
        elem_shock_dir = int(input('Provide element count along shock direction: '))
        elem_perp_1 = int(input('Provide element count perpendicular to shock direction: '))
        elem_perp_2 = 1

        #format and replace

        final = final.replace("{ELEM_PERP_1}", str(elem_perp_1))
        final = final.replace("{ELEM_PERP_2}", str(elem_perp_2))
        final = final.replace("{ELEM_SHOCK_DIR}", str(elem_shock_dir))

        final = final.replace("{DIM_PERP_1}", str(float(elemsize[-1] * elem_perp_1)))
        final = final.replace("{DIM_PERP_2}", str(float(elemsize[-1] * elem_perp_2) if pbx != 'LOADED' else str(zloc + elemsize[-1])))
        final = final.replace("{DIM_SHOCK_DIR}", str(float(elemsize[0] * elem_shock_dir)))
        
        #for the case wher we want a slice from a loaded 3D structure
        if (pbx == 'LOADED' and dim == '2D'):
            repl = str(zloc)
        else:
            repl = str(0)

        final = final.replace("{ZMIN}", repl)

    case '3D':
        elem_shock_dir = int(input('Provide element count along shock direction: '))
        elem_perp_1 = int(input('Provide element count perpendicular to shock direction: '))
        elem_perp_2 = elem_perp_1 #square cross section

        final = final.replace("{ELEM_PERP_1}", str(elem_perp_1))
        final = final.replace("{ELEM_PERP_2}", str(elem_perp_2))
        final = final.replace("{ELEM_SHOCK_DIR}", str(elem_shock_dir))

        final = final.replace("{DIM_PERP_1}", str(float(elemsize[-1] * elem_perp_1)))
        final = final.replace("{DIM_PERP_2}", str(float(elemsize[-1] * elem_perp_2)))
        final = final.replace("{DIM_SHOCK_DIR}", str(float(elemsize[0] * elem_shock_dir)))

        #for the case wher we want a slice from a loaded 3D structure
        if (pbx == 'LOADED' and dim == '2D'):
            repl = str(zloc)
        else:
            repl = str(0)
            
        final = final.replace("{ZMIN}", repl)

#define name for the loaded microstructures
if (pbx == 'LOADED'):
    final_name = f"loaded_{datafile}_slice{zloc}_{dim}_up{vel}_perp{elem_perp_1}_shock{elem_shock_dir}_time{'distr' if use_distributions == 'true' else 'tabular'}_binder{binder_width}_{int(binder_range[0])}_{int(binder_range[-1])}" 
else:
    final_name = f"{dim}_up{vel}_type{pbx}_perp{elem_perp_1}_shock{elem_shock_dir}_poro{particle_poro if pbx == 'PBX' else 0}_time{'distr' if use_distributions == 'true' else 'tabular'}_binder{binder_width}_{int(binder_range[0])}_{int(binder_range[-1])}" 

with open(f'{final_name}.i', "w") as f:
    f.write(final)

print('File Generated !!')
print('#################')

gen_and_run = input('Generate sbatch script and submit? (YES) or (NO): ')

if (gen_and_run == 'YES'):
    nodes, cores, days, hours, app = input('Provide nodes, cores, days, hours, and app name: ').split(',')
    #generate a new dir for every run
    su.run(['mkdir', f'DIR{final_name}'])

    #copy all csv files to each particular folder
    import glob
    csv_files = glob.glob('csv/*')
    su.run(['cp', *csv_files, f'DIR{final_name}'], check=True)

    #move the input to that directory
    su.run(['mv', f'{final_name}.i', f'DIR{final_name}'])

    #go to the created directory
    su.run(['cd', f'DIR{final_name}'])

    #create sbatch file inside the created directory
    su.run(['touch', f'DIR{final_name}/{final_name}'])

    #paste the .txt file if required
    if (pbx == 'LOADED'):
        su.run(['cp', f'{datafile}.txt', f'DIR{final_name}'])

    #use sbatch template to generate actual sbatch file within this directory
    with open(f'sbatch_template', 'r') as f:
        sbatch_template = f.read()

    sbatch_final = sbatch_template.replace('{NODES}', str(nodes))
    sbatch_final = sbatch_final.replace('{CPUS}', str(cores))
    sbatch_final = sbatch_final.replace('{DAYS}', str(days))
    sbatch_final = sbatch_final.replace('{HOURS}', str(hours))
    sbatch_final = sbatch_final.replace('{APP_NAME}', str(app))
    sbatch_final = sbatch_final.replace('{INPUT}', f'{final_name}')
    sbatch_final = sbatch_final.replace('{DIR}', f'DIR{final_name}')

    #input queue
    queue = str(input('Queue to submit the job (normal) or (standby): '))
    sbatch_final = sbatch_final.replace('{QUEUE}', str(queue))

    with open(f'DIR{final_name}/{final_name}', 'w') as f:
        f.write(sbatch_final)
    
    #submit the job
    su.run(['sbatch', f'DIR{final_name}/{final_name}'])

#generate log file with al the used configuration
dump_name = f'DIR{final_name}/DUMP{final_name}.txt'
with open(dump_name, 'w') as f:
    f.write(f"Simulation Dimension: {dim}\n")
    if dim == '2D':
        f.write(f"Element Count along Shock Direction: {elem_shock_dir}\n")
        f.write(f"Element Count Perpendicular to Shock Direction: {elem_perp_1}\n")
    else:
        f.write(f"Element Count along Shock Direction: {elem_shock_dir}\n")
        f.write(f"Element Count Perpendicular to Shock Direction: {elem_perp_1} (square cross section)\n") 
    f.write(f"Impact Velocity: {vel} km/s\n")
    f.write(f"Microstructure Type: {pbx}\n")
    if pbx == 'PBX':
        f.write(f"Particle Sizes (microns): {particle_sizes}\n")
        f.write(f"Particle Proportions: {particle_proportions}\n")
        f.write(f"Particle Porosity: {particle_poro*100}%\n")
        f.write(f"Binder Width (microns): {binder_width}\n")
        f.write(f"Range of Binder Microstructure IDs: {binder_range}\n")
    else:
        f.write(f"Range of Microstructure IDs: {range_microstructures}\n")

    f.write(f"Bulk Grains?: {bulk_grains}\n")
    f.write(f"Use Mixture for Mechanics?: {use_mixture}\n")
    f.write(f"Binder Properties (Bulk, Yield, Shear Modulus, Poisson's Ratio): {binder_properties}\n")
    f.write(f"Use Tabular Time?: {tabular_time}\n")
    f.write(f"Use Time Distributions?: {use_distributions}\n")
    f.write(f"Use Gating for Heat and Chemical Sources?: {use_gating}\n")
#################################################################################
