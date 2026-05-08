import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#settings for plots

plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
})

#helper functions for error estimation
us = lambda up: -0.286 * up**3 + 1.640 * up**2 - 1.249391 * up**1 + 5.575975

#define function to estimate error
def estimError(up_avg, h):
    up_avg = up_avg / 100 #scaling factor
    error = 0
    #intial shock velocity
    us_i = us(up_avg)

    #1. Element size
    error += h

    #2. Shock width at up during detonation
    error += 5 * h

    #3.  Error estimation in timestep resolution
    dt = 5e-3 
    period = 20
    dt_real = period * dt

    #distance traveled by shock at the given resolution
    dist = us_i * dt_real
    error += dist
    return error

#helper function for bisection
def find_intersection(fpre, fpost, maxiter, tol):
    #define the target function for the bisection
    target = lambda x: fpre(x) - fpost(x)
    #run a bisection to find the roots
    #initial function evaluations
    l = 0
    r = 1e3
    m = (l + r) / 2
    fl = target(l)
    fr = target(r)
    fm = target(m)
    it = 0
    while it < maxiter:
        it += 1
        if (fl * fm) < 0:
            l = l
            r = m
            m = (l + r) / 2
            fl = fl
            fr = fm
            fm = target(m)
        else:
            l = m
            r = r
            m = (l + r) / 2
            fl = fm
            fr = fr
            fm = target(m)
        #check convergence
        if (abs(fm) < tol):
            print(f'Convergence achieved. Detonation distance {m}')
            return m
    print(f'Diverged. Last value {m}')
    return m

datasets = [175,185,195]

#import the new csv data
data_dict = {}

def find_run_det(data):
    #for a position in y, time in x axis plot
    xdata = data['Points:1']
    ydata = data['tracking']

    #flip y data
    ydata = np.abs(ydata - ydata.max())
    
    #smooth the data
    ydata = np.convolve(ydata, np.ones(20)/20, mode='same')

    #compute gradient
    grad = np.gradient(ydata, xdata)

    #take derivative of gradient
    grad2 = np.gradient(grad, xdata)

    #form mask for gradient
    mask = (xdata < 250) & (xdata > 20)
    
    #find biggest change in grad
    interval = np.where(mask)[0]
    max_grad_index = interval[np.argmax(np.abs(grad2[mask]))]

    #compute the average value of the gradient before and after detonation
    #this is to get the two lines for the pre and post detonation

    avg_grad_pre = np.mean(grad[(xdata < xdata[max_grad_index]) & mask])
    avg_grad_post = np.mean(grad[(xdata > xdata[max_grad_index]) & mask])

    #form linear functions with these gradients
    fpre_avg = lambda x: avg_grad_pre * (x - xdata[max_grad_index]) + ydata[max_grad_index]
    fpost_avg = lambda x: avg_grad_post * (x - xdata[max_grad_index]) + ydata[max_grad_index]

    #take every pair of points and get the slope before and after
    #this will generate a distribution of slopes, we can take the average and standard deviation to get an estimate of the error in the position of detonation
    slopes_pre = []
    slopes_post = []
    for i in range(len(xdata)):
        for j in range(i+1, len(xdata)):
            if (xdata[i] < xdata[max_grad_index]) and (xdata[j] < xdata[max_grad_index]) and mask[i] and mask[j]:
                slope = (ydata[j] - ydata[i]) / (xdata[j] - xdata[i])
                slopes_pre.append(slope)

    for i in range(len(xdata)):
        for j in range(i+1, len(xdata)):
            if (xdata[i] > xdata[max_grad_index]) and (xdata[j] > xdata[max_grad_index]) and mask[i] and mask[j]:
                slope = (ydata[j] - ydata[i]) / (xdata[j] - xdata[i])
                slopes_post.append(slope)

    #obtain the averages of the dsitribution
    avg_slope_pre = np.mean(slopes_pre)
    avg_slope_post = np.mean(slopes_post)

    #obtain the standard deviation of the distribution
    std_slope_pre = np.std(slopes_pre)
    std_slope_post = np.std(slopes_post)

    #maximum and minimum
    max_slope_pre = np.max(slopes_pre)
    min_slope_pre = np.min(slopes_pre)

    max_slope_post = np.max(slopes_post)
    min_slope_post = np.min(slopes_post)

    #form maximum and minimum lines
    fpre_max = lambda x: max_slope_pre * (x - xdata[max_grad_index]) + ydata[max_grad_index]    
    fpre_min = lambda x: min_slope_pre * (x - xdata[max_grad_index]) + ydata[max_grad_index]

    fpost_max = lambda x: max_slope_post * (x - xdata[max_grad_index]) + ydata[max_grad_index]
    fpost_min = lambda x: min_slope_post * (x - xdata[max_grad_index]) + ydata[max_grad_index]

    #form line for standard deviation case
    fpre_std_max = lambda x: (avg_slope_pre + std_slope_pre) * (x - xdata[max_grad_index]) + ydata[max_grad_index]
    fpre_std_min = lambda x: (avg_slope_pre - std_slope_pre) * (x - xdata[max_grad_index]) + ydata[max_grad_index]

    #find the corresponding position and time
    run_det = xdata[max_grad_index]
    time_det = ydata[max_grad_index]
    return run_det, time_det, xdata, ydata, fpre_avg, fpost_avg, avg_slope_pre, std_slope_pre, avg_slope_post, std_slope_post, fpre_max, fpre_min, fpost_max, fpost_min, fpre_std_max, fpre_std_min

rmin, rmax, ravg, rstd, times = [], [], [], [], []

#append other data for the pop plot

fig, ax = plt.subplots()
figpop, axpop = plt.subplots()
figstat, axstat = plt.subplots(figsize = (10,7))

#simulation run data
Ps_sim = [19, 23, 28]

for i, dataset in enumerate(datasets):
    data_dict[dataset] = pd.read_csv(f'{dataset}_tracking.csv')
    run_det, time_det, xdata, ydata, fpre, fpost, avg_slope_pre, std_slope_pre, avg_slope_post, std_slope_post, fpre_max, fpre_min, fpost_max, fpost_min, fpre_std_max, fpre_std_min = find_run_det(data_dict[dataset])
    times.append(time_det)

    #compute intersection of max, min, and average lines
    intersection_avg = find_intersection(fpre, fpost, 1e4, 1e-5)
    intersection_max = find_intersection(fpre_max, fpost_min, 1e4, 1e-5)
    intersection_min = find_intersection(fpre_min, fpost_max, 1e4, 1e-5)
    intersection_std_max = find_intersection(fpre_std_max, fpost, 1e4, 1e-5)
    intersection_std_min = find_intersection(fpre_std_min, fpost, 1e4, 1e-5)

    rmin.append(intersection_min)
    rmax.append(intersection_max)
    ravg.append(intersection_avg)
    rstd.append(abs(intersection_std_max - intersection_std_min))

    #plot distributions including max, min, avg, and standard deviation
    axstat.errorbar(intersection_avg, Ps_sim[i], xerr = estimError(float(dataset)/100, 1.9) + np.abs(intersection_avg - intersection_min), fmt='ko', label=f'{dataset}')
    axstat.set_xlabel('Run to Detonation $(\mu m)$', fontsize=20)
    axstat.set_ylabel('Impact Pressure (GPa)', fontsize=20)
    axstat.tick_params(axis='both', labelsize=20)

    #plot 
    ax.plot(xdata, ydata, label=f'{float(dataset)/100}, detonation onset {rmin[i]:.2f} $\\mu m$')
    ax.plot(run_det, time_det, 'ko')
    ax.set_xlabel('Run Distance $(\mu m)$', fontsize=20)
    ax.set_ylabel('Run Time (ns)', fontsize=20)
    ax.tick_params(axis='both', labelsize=20)
    ax.plot(xdata, fpre(xdata), 'k--')
    ax.plot(xdata, fpost(xdata), 'k--')
    ax.legend(fontsize=20)

    ax.set_xlim(30, 120)
    ax.set_ylim(2, 20)
    ax.grid()

    axpop.plot(run_det, time_det, 'o', label=f'{ravg[i]} $\\mu m$')
    axpop.set_xlabel('Run Distance $(\mu m)$', fontsize=20)
    axpop.set_ylabel('Run Time (ns)', fontsize=20)
    axpop.tick_params(axis='both', labelsize=20)
    axpop.legend(fontsize=20)
    axpop.set_xlim(1e1, 1e3)

plt.show()

#now form the the pop plot with the errorbar data
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy.optimize as op

#define functions for each material

pop_9407 = lambda P: np.exp((0.57 - np.log(P)) / 0.49)
range_9407 = np.linspace(1.14, 4.69, 100) #GPa

#ELVAX
pop_elvax = lambda P: np.exp((1.43 - np.log(P)) / 0.73)
range_elvax = np.linspace(2.96, 11.71, 100) #GPa

#pop_CompB = lambda P: 1e3 * 2000 * P ** (-1.34) #pressure in kbar
#range_CompB = np.linspace(4, 10, 100) #GPa

##jason fittings
pop_jason1 = lambda P: np.exp(1.71 - 1.43 * np.log(P))
range_jason1 = range_elvax

pop_HMX = lambda P: np.exp((1.18 - np.log(P)) / 0.59)
range_HMX = np.linspace(4.41, 9.55, 100) #GPa

pop_9501 = lambda P: np.exp(1.15 - np.log(P) / 0.64)
range_9501 = np.linspace(2.38, 7.32, 100)

pop_LX4 = lambda P: np.exp(1.01 - np.log(P) / 0.47)
range_LX4 = np.linspace(4.06, 6.74, 100)

#the runs are used from the initial part of the code
x_sim = [value for value in ravg]
x_sim_err = [value for value in rstd]

#power law fit function for data to look linear in log-log plot

def fitfunc(P, A, n):
    return A * P ** n

#do function fitting
from scipy.optimize import curve_fit
fitted_params = curve_fit(fitfunc, x_sim, Ps_sim)

#take the fitted parameters and form the function for evaluation
Afit, nfit = fitted_params[0]
zsim = lambda P: Afit * P ** nfit

plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
})

fig, ax = plt.subplots(figsize=(15,13))
import matplotlib.cm as cm
import matplotlib.colors as mcolors

cmap = cm.get_cmap('jet')
norm = mcolors.Normalize(vmin=0, vmax=1)
cmap_norm = cm.ScalarMappable(norm=norm, cmap=cmap)

ax.loglog(pop_9407(range_9407), range_9407, c = cmap(norm(0.1)), label = 'PBX 9407 (94% RDX, 6% Exon)', linewidth = 5)
#ax.loglog(pop_CompB(range_CompB), range_CompB, c = cmap(norm(0.4)), label = 'CompB (65% RDX, 35% TNT)')
ax.loglog(pop_HMX(range_HMX), range_HMX, c = cmap(norm(0.3)), label = 'HMX (100% HMX)', linewidth = 5)
ax.loglog(pop_9501(range_9501), range_9501, c = cmap(norm(0.5)), label = 'PBX 9501 (95% HMX, 5% Estane)', linewidth = 5)
#ax.loglog(pop_LX4(range_LX4), range_LX4, c = cmap(norm(0.8)), label = 'LX4 (85% HMX, 15% Viton)')
ax.loglog(pop_elvax(range_elvax), range_elvax, c = cmap(norm(0.7)), label = 'RDX/2.5 WAX/2.5 ELVAX \n (95% RDX, 2.5% WAX/ELVAX)', linewidth = 5)
##ax.errorbar(x_sim[0], Ps_sim[0],(22.18 - 21.11),(95.38-90.11),fmt='k*', label = 'Simulation Data')
ax.set_xlabel('Run to Detonation (mm)', fontsize = 50)
ax.set_ylabel('Initial Impact Pressure (GPa)', fontsize = 50)
ax.loglog(1e-3 * np.linspace(x_sim[0], x_sim[-1], 100), zsim(np.linspace(x_sim[0], x_sim[-1], 100)), 'k--', linewidth = 5)
ax.loglog(pop_jason1(range_jason1), range_jason1, c = cmap(norm(0.9)), label = 'PBX 9501 (95% HMX, 5% Estane)', linewidth = 5)

#box = ax.get_position()
#ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

#plot the error bars
for i in range(len(x_sim)):
    ax.loglog(1e-3 * x_sim[i], Ps_sim[i], 'k*', markersize = 20, label = 'Simulation Data' if i == 0 else "")

plt.xticks(fontsize = 50)
plt.yticks(fontsize = 50)
#ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), ncol = 2, fontsize=20)
#ax.legend(loc='lower left', ncol = 1, fontsize=30, frameon=False)
for spine in ax.spines.values():
    spine.set_linewidth(2)
#plt.tight_layout()
plt.savefig('updated_pop.jpeg', dpi=500)
plt.show()