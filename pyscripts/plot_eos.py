import numpy as np
from paraview.simple import *
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from vtk.util import numpy_support
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pandas as pd

#make
#retrieve inputs
filename = input('Provide the file name without extension: ')
num_slices = int(input('Provide number of slices in space for sampling: '))
min_time = int(input('Provide minimum timestep tag for analysis: '))
max_time = int(input('Provide maximun timestep tag for analysis: '))
num_times = int(input('Provide number of time points to analyze: '))
xaxis_var = input('x axis EOS variable: ')
xaxis_type = input('x axis type: ')
yaxis_var = input('y axis EOS variable: ')
yaxis_type = input('y axis type: ')
plot_fit = input('Plot fitted curve? (YES) or (NO) ')
print_distributions = input('Print distributions? ')
color_by = input('Color data by: ').strip()
color_mode = color_by.upper()

if (plot_fit == 'YES'):
    omega_unreacted = float(input('omega value for the unreacted branch: '))
    omega_reacted = float(input('omega value for the reacted branch: '))
    base_unreacted = f"DATA_unreacted_eos_BOTH_omegaU{omega_unreacted}_omegaR{omega_reacted}.csv"
    base_reacted = f"DATA_reacted_eos_BOTH_omegaU{omega_unreacted}_omegaR{omega_reacted}.csv"

    #read data from csv files
    df_unreacted = pd.read_csv(base_unreacted, sep=',')
    df_reacted = pd.read_csv(base_reacted, sep=',')

#allocate slice locations
nslices = num_slices
size = 380
locs = [(size / nslices)*t for t in range(nslices)]

#generate location vectors
positions = []
for loc in locs:
    positions.append([27.5, loc, 0.55])

#form file reader
fname = f'../{filename}_out.e'
reader = ExodusIIReader(FileName=[fname])
reader.UpdatePipeline()

#generate slicer and direction
slicer = Slice(Input=reader)
slicer.SliceType = 'Plane'
slicer.SliceType.Normal = [0, 1, 0]

#define times
dt = num_times
ts = np.asarray(reader.TimestepValues, dtype = float)
idx = np.linspace(min_time, max_time, dt).astype(int)
times = ts[idx]

#define arrays
mins = np.full((len(times), len(positions)), np.nan)
maxs = np.full((len(times), len(positions)), np.nan)
Q1 = np.full((len(times), len(positions)), np.nan)
Q3 = np.full((len(times), len(positions)), np.nan)
avg = np.full((len(times), len(positions)), np.nan)
std = np.full((len(times), len(positions)), np.nan)

xavg = np.full((len(times), len(positions)), np.nan)

#for coloring
coloravg = np.full((len(times), len(positions)), np.nan)
colormin = np.full((len(times), len(positions)), np.nan)
colormax = np.full((len(times), len(positions)), np.nan)

def split_composite(array):
    if hasattr(array, "Arrays"):
        parts = []
        for blockarr in array.Arrays:
            if blockarr is None:
                continue
            parts.append(np.asarray(blockarr, dtype=float).ravel())
        res = np.concatenate(parts) if parts else np.array([])
    else:
        res = np.asarray(array, dtype=float).ravel()
        
    res = res[np.isfinite(res)]
    return res

for it, t in enumerate(times): 
    for ix, pos in enumerate(positions):
        slicer.SliceType.Origin = pos
        slicer.UpdatePipeline(time=t)
        
        #compute statistics
        data_raw = dsa.WrapDataObject(servermanager.Fetch(slicer))
        if (xaxis_type == 'POINT' or yaxis_type == 'POINT'):
            point_raw = data_raw.PointData
            
        cell_raw = data_raw.CellData
        
        #get point and cell data
        if (xaxis_type == 'POINT'):
            xdata = point_raw[xaxis_var]
        elif (xaxis_type == 'CELL'):
            xdata = cell_raw[xaxis_var]
            
        if (yaxis_type == 'POINT'):
            ydata = point_raw[yaxis_var]
            #colordata = point_raw[color_by]

        elif (yaxis_type == 'CELL'):
            ydata = cell_raw[yaxis_var]
            #colordata = cell_raw[color_by]
        
        #split array
        x_split = split_composite(xdata)
        y_split = split_composite(ydata)
        #color_split = split_composite(colordata)

        if color_mode == 'TIME':
            color_split = None
        else:
            if(yaxis_type == 'POINT'):
                color_data = point_raw[color_by]
            elif(yaxis_type == 'CELL'):
                color_data = cell_raw[color_by]
            color_split = split_composite(color_data)

        if y_split.size == 0:
            continue
        else:
            mins[it, ix] = float(y_split.min())
            maxs[it, ix] = float(y_split.max())
            Q1[it, ix] = float(np.percentile(y_split, 25))
            Q3[it, ix] = float(np.percentile(y_split, 75))
            avg[it, ix] = float(y_split.mean())
            std[it, ix] = float(y_split.std())

            if color_mode == 'TIME':
                coloravg[it, ix] = float(t)
                colormin[it, ix] = float(t)
                colormax[it, ix] = float(t)
            else:
                coloravg[it, ix] = float(color_split.mean())
                colormin[it, ix] = float(color_split.min())
                colormax[it, ix] = float(color_split.max())
                
        if x_split.size == 0:
            continue
        else:
            xavg[it, ix] = float(x_split.mean())
        

######PLOT######
import matplotlib.cm as cm
import matplotlib.colors as mcolors

plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
})

fig, ax = plt.subplots()
cmap = cm.get_cmap('jet')
#norm = mcolors.Normalize(vmin=float(colormin.min()), vmax=float(colormax.max()))
#sm = cm.ScalarMappable(norm=norm, cmap=cmap)

#for coloring

if color_mode == 'TIME':
    norm = mcolors.Normalize(vmin=float(times.min()), vmax=float(times.max()))
else:
    all_colors = coloravg[np.isfinite(coloravg)]
    norm = mcolors.Normalize(vmin=float(all_colors.min()), vmax=float(all_colors.max()))

sm = cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

#pressures
ps = np.array([p[0] for p in positions], dtype = float)
for it, t in enumerate(times):
    
    validavg = np.isfinite(xavg[it, :]) & np.isfinite(avg[it, :]) & np.isfinite(coloravg[it, :])
    validmin = np.isfinite(xavg[it, :]) & np.isfinite(avg[it, :]) & np.isfinite(colormin[it, :])
    validmax = np.isfinite(xavg[it, :]) & np.isfinite(avg[it, :]) & np.isfinite(colormax[it, :])
    
    ax.scatter(
        xavg[it, validavg],
        np.abs(avg[it, validavg]),
        c=coloravg[it, validavg],
        linewidth=2,
        marker='o',
        cmap=cmap,
        norm=norm
    )

    if (print_distributions == 'YES'):
        #ax.scatter(xavg[it, :], (Q1[it, :]), color=color, linewidth = 1, alpha=0.5, marker = '<')
        #ax.scatter(xavg[it, :], (Q3[it, :]), color=color, linewidth = 1, alpha=0.5, marker = '>')
        ax.scatter(xavg[it, validmin], (mins[it, validmin]), c=colormin[it, validmin], linewidth = 0.5, alpha=0.25, marker = 'v', cmap=cmap)
        ax.scatter(xavg[it, validmax], (maxs[it, validmax]), c=colormax[it, validmax], linewidth = 0.5, alpha=0.25, marker = '^', cmap=cmap)
        #ax.plot(xavg[it, :], (std[it, :]), color=color, linewidth = 0.75, alpha=0.3, marker='+')

#plot fitted curve
if (plot_fit == 'YES'):
    ax.plot(df_unreacted['J_unreacted'].to_numpy(), df_unreacted['P_unreacted'], 'k-', linewidth = 1)
    ax.plot(df_reacted['J_reacted'].to_numpy(), df_reacted['P_reacted'], 'k-', linewidth = 1)

##for legend only
if (print_distributions == 'YES'):
    ax.plot([],[], 'ko', label='Average values')
    #ax.plot([],[], 'k<', label='Q1 values')
    #ax.plot([],[], 'k>', label='Q3 values')
    ax.plot([],[], 'kv', label='Minimum values')
    ax.plot([],[], 'k^', label='Maximum values')
    #ax.plot([],[], '+', label='Standard deviation')

ax.set_xlabel(f'{xaxis_var}', fontsize=20)
ax.set_ylabel(f'{yaxis_var}', fontsize=20)
#ax.set_xlim(1.5, 3.6)
#ax.set_ylim(0, 60)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

cbar = fig.colorbar(sm, ax=ax)
cbar.set_label(f'{color_by}', fontsize=20)
plt.tight_layout()
#plt.grid(True)
plt.legend()
plt.savefig(f'{xaxis_var}{yaxis_var}_eos_{num_slices}_samples.jpeg', dpi=1000, bbox_inches='tight')
plt.show()
