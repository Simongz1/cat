import numpy as np
from paraview.simple import *
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from vtk.util import numpy_support
import matplotlib.pyplot as plt

#retrieve inputs
filename = input('Provide the file name without extension: ')
dimension = input('Provide the dimension of the simulation (2D) or (3D): ')
size = float(input('Domain size along the shock direction: '))
trans = float(input('Domain size perpendicular to the shock direction: '))
field = input('Field to plot: ')
vartype = input('Field type for elemenal (CELL) or nodal (POINT) association: ')
num_slices = int(input('Provide number of slices in space: '))
min_time = int(input('Provide minimum timestep tag for analysis: '))
max_time = int(input('Provide maximun timestep tag for analysis: '))
num_times = int(input('Provide number of time history lines to generate: '))
include_stats = input('Include statistics (YES) or (NO): ')

#allocate slice locations
nslices = num_slices
locs = np.linspace(0, size, nslices + 1)

#generate location vectors
#positions = []
#for loc in locs:
#    positions.append([loc, 0, 0])

#form file reader
fname = f'../{filename}_out.e'
reader = ExodusIIReader(FileName=[fname])
reader.UpdatePipeline()

#generate slicer and direction
match dimension:
    case '3D':
        slicer = Slice(Input=reader)
        slicer.SliceType = 'Plane'
        slicer.SliceType.Normal = [0, 1, 0]
        positions = []
        for loc in locs:
            positions.append([27.5, loc, 0.55])
    case '2D':
        slicer = Clip(Input=reader)
        slicer.ClipType = 'Box'
        slicer.ClipType.Length = [trans, 10, 1.1]
        positions = []
        for loc in locs:
            positions.append([0, loc, 0])

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

def split_composite(array):

    if array is None:
        return np.array([])
    if array.__class__.__name__ == "VTKNoneArray":
        return np.array([])
    
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
        match dimension:
            case '2D':
                slicer.ClipType.Position = pos
                slicer.UpdatePipeline(time=t)

            case '3D':
                slicer.SliceType.Origin = pos
                slicer.UpdatePipeline(time=t)
        
        #compute statistics
        data_raw = dsa.WrapDataObject(servermanager.Fetch(slicer))
        if (vartype == 'POINT'):
            data_raw = data_raw.PointData
        elif (vartype == 'CELL'):
            data_raw = data_raw.CellData

        data_raw = data_raw[field]
        
        #split array
        data_split = split_composite(data_raw)
        if data_split.size == 0:
            match dimension:
                case '2D':
                    mins[it, ix] = np.nan
                    maxs[it, ix] = np.nan
                    Q1[it, ix]   = np.nan
                    Q3[it, ix]   = np.nan
                    avg[it, ix]  = np.nan
            continue
        else:
            mins[it, ix] = float(data_split.min())
            maxs[it, ix] = float(data_split.max())
            Q1[it, ix] = float(np.percentile(data_split, 25))
            Q3[it, ix] = float(np.percentile(data_split, 75))
            avg[it, ix] = float(data_split.mean())
            
            
######PLOT######
import matplotlib.cm as cm
import matplotlib.colors as mcolors

plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
})

fig, ax = plt.subplots()
cmap = cm.get_cmap('jet')
norm = mcolors.Normalize(vmin=float(times.min()), vmax=float(times.max()))
sm = cm.ScalarMappable(norm=norm, cmap=cmap)

#pressures
ps = np.array([p[1] for p in positions], dtype = float)
for it, t in enumerate(times):
    color = cmap(norm(float(t)))
    
    ax.plot(ps, np.abs(avg[it, :]), color='k', linewidth = 2)
    if (include_stats == 'YES'):
        ax.plot(ps, np.abs(Q1[it, :]), color='k', linewidth = 0.5, alpha=0.35)
        ax.plot(ps, np.abs(Q3[it, :]), color='k', linewidth = 0.5, alpha=0.35)
        ax.plot(ps, np.abs(mins[it, :]), color='k', linewidth = 0.5, alpha=0.35)
        ax.plot(ps, np.abs(maxs[it, :]), color='k', linewidth = 0.5, alpha=0.35)
        
        ax.fill_between(ps, np.abs(Q1[it, :]), np.abs(Q3[it, :]), color = color, alpha=0.25, linewidth=0)
        ax.fill_between(ps, np.abs(mins[it, :]), np.abs(maxs[it, :]), color = color, alpha=0.1, linewidth=0)
        
ax.set_xlabel('Distance($\mu m$)', fontsize=20)
ax.set_ylabel(f'{field}', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

cbar = fig.colorbar(sm, ax=ax)
cbar.set_label('time (ns)', fontsize=20)
plt.tight_layout()
plt.savefig(f'wide{dimension}_{field}_stats_{num_slices}slices_{min_time}-{max_time}.jpeg', dpi=1000, bbox_inches='tight')
plt.show()
