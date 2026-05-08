import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy.optimize as op

#user inputs
branch = input('(REACTED), (UNREACTED), (BOTH)?: ')

#request name of images for T-V and P-V data
TV = input('Absolute path for T-V image: ')
PV = input('Absolute path for P-V image: ')
match branch:
    case "UNREACTED":
        omega = float(input('Estimation for unreacted omega: '))
    case "REACTED":
        omega = float(input('Estimation for reacted omega: '))
    case "BOTH":
        omega_unreacted = float(input('Estimation for unreacted omega: '))
        omega_reacted = float(input('Estimation for reacted omega: '))
        

def ScaleArray(array, xmin, xmax, ymin, ymax):
    x = np.asarray(array, dtype=float)
    new = (x - xmin) / (xmax - xmin) * (ymax - ymin) + ymin
    return new

#extract data from the image using ginput

def Normalize(image, min_x, max_x, min_y, max_y):
    #get min and max for x and y axis:
    image = mpimg.imread(image)
    show = plt.imshow(image)
    (xmin, ymin), = plt.ginput(1)
    (xmax, ymax), = plt.ginput(1)

    #points
    points = plt.ginput(0, 0)

    plt.close()

    #data
    x = np.array([p[0] for p in points])
    y = np.array([p[1] for p in points])

    #ranges:
    maxx, minx = max_x, min_x
    maxy, miny = max_y, min_y

    #scale data
    x_scaled = ScaleArray(x, xmin, xmax, minx, maxx)
    y_scaled = ScaleArray(y, ymin, ymax, miny, maxy)

    return x_scaled, y_scaled

match branch:
    case "UNREACTED" | "REACTED":  
        JT_scaled, T_scaled = Normalize('RDX_TV.png', 0.6, 1, 0, 4000)

match branch:
    case "BOTH":
        JT_scaled_unreacted, T_scaled_unreacted = Normalize(TV, 0.6, 1, 0, 4000)
        JT_scaled_reacted, T_scaled_reacted = Normalize(TV, 0.6, 1, 0, 4000)

#label name for force field used
ff = input('Force field name: ')

#compute the Pth values
def getPth(omega, J, T):
    rho = 1.820 #originally was 1.5
    cv = 2320e-6
    values = [omega * rho * cv * T[i] / J[i] for i in range(len(T))]
    return values

#get thermal pressure part
match branch:
    case "UNREACTED" | "REACTED":  
        Ph = getPth(omega, JT_scaled, T_scaled)

match branch:
    case "BOTH":
        Ph_unreacted = getPth(omega_unreacted, JT_scaled_unreacted, T_scaled_unreacted)
        Ph_reacted = getPth(omega_reacted, JT_scaled_reacted, T_scaled_reacted)
        

#polynomial fit of 3rd degree for Ph - J relation
match branch:
    case "UNREACTED":
        degree = 4
        fitPh = np.polyfit(JT_scaled, Ph, degree)
        fitPh = np.poly1d(fitPh) #this function gets J-Ph relation as a callable
    case "REACTED":
        degree = 1
        fitPh = np.polyfit(JT_scaled, Ph, degree)
        fitPh = np.poly1d(fitPh) #this function gets J-Ph relation as a callable
    case "BOTH":
        degree_unreacted = 4
        degree_reacted = 1

        fitPh_unreacted = np.polyfit(JT_scaled_unreacted, Ph_unreacted, degree_unreacted)
        fitPh_reacted = np.polyfit(JT_scaled_reacted, Ph_reacted, degree_reacted)

        fitPh_unreacted = np.poly1d(fitPh_unreacted) #this function gets J-Ph relation as a callable
        fitPh_reacted = np.poly1d(fitPh_reacted) #this function gets J-Ph relation as a callable

#get actual Pv data
match branch:
    case "UNREACTED" | "REACTED":  
        JP_scaled, P_scaled = Normalize('RDX_PV.png', 0.6, 1, 0, 60)

match branch:
    case "BOTH":
        JP_scaled_unreacted, P_scaled_unreacted = Normalize(PV, 0.6, 1, 0, 60)
        JP_scaled_reacted, P_scaled_reacted = Normalize(PV, 0.6, 1, 0, 60)
        

#compute cold pressure
match branch:
    case "UNREACTED" | "REACTED":
        Pc = [P_scaled[i] - fitPh(JP_scaled[i]) for i in range(len(JP_scaled))]

match branch:
    case "BOTH":
        Pc_unreacted = [P_scaled_unreacted[i] - fitPh_unreacted(JP_scaled_unreacted[i]) for i in range(len(JP_scaled_unreacted))]
        Pc_reacted = [P_scaled_reacted[i] - fitPh_reacted(JP_scaled_reacted[i]) for i in range(len(JP_scaled_reacted))]
        

#fit the cold curve agains the cold part of JWL

def coldJWL(J, A, B, R1):
    return A * np.exp(-R1 * J) + B * np.exp(- 0.1 * R1 * J)

match branch:
    case "UNREACTED":
        p0 = [2.3e3, -1e2, 1e1]
    case "REACTED":
        p0 = [1e3, 1e3, 1e2]

    case "BOTH":
        p0_unreacted = [2.3e3, -1e2, 1e1]
        p0_reacted = [1e3, 1e3, 1e2]

match branch:
    case "UNREACTED" | "REACTED":
        popt, pcov = op.curve_fit(coldJWL, JP_scaled, Pc, p0 = p0, maxfev = 10000)

match branch:
    case "BOTH":
        popt_unreacted, pcov_unreacted = op.curve_fit(coldJWL, JP_scaled_unreacted, Pc_unreacted, p0 = p0_unreacted, maxfev = 10000)
        popt_reacted, pcov_reacted = op.curve_fit(coldJWL, JP_scaled_reacted, Pc_reacted, p0 = p0_reacted, maxfev = 10000)

#reconstruct total pressure
fig, ax = plt.subplots()

Jeval = np.linspace(0.6, 1, 100)
plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
})

match branch:
    case "UNREACTED" | "REACTED":
        P_total = [coldJWL(J, *popt) + fitPh(J) for J in Jeval]
        ax.plot(Jeval, P_total, 'k-')
        ax.plot(JP_scaled_unreacted, P_scaled_unreacted, 'r*')
        ax.set_title(f'JWL Paramters for the {branch} branch are: A = {popt[0]:.4f} \n B = {popt[1]:.4f}, R1 = {popt[2]:.4f} \n R2 = {0.1 * popt[2]:.4f}')

        plt.savefig(f'eosimg/eos_{branch}_A{popt[0]:.4f}_B{popt[1]:.4f}_R1{popt[2]:.4f}_R2{0.1 * popt[2]:.4f}_omega{omega}.jpeg', dpi = 1500, bbox_inches='tight')

match branch:
    case "BOTH":
        P_total_unreacted = [coldJWL(J, *popt_unreacted) + fitPh_unreacted(J) for J in Jeval]
        P_total_reacted = [coldJWL(J, *popt_reacted) + fitPh_reacted(J) for J in Jeval]

        ax.plot(Jeval, P_total_unreacted, 'b-')
        ax.plot(Jeval, P_total_reacted, 'r-')

        ax.plot(JP_scaled_unreacted, P_scaled_unreacted, 'k*')
        ax.plot(JP_scaled_reacted, P_scaled_reacted, 'k*')

        base_name = f'FF{ff}_eos_{branch}_omegaU{omega_unreacted}_omegaR{omega_reacted}'
        #set tick sizes
        ax.tick_params(axis='both', which='major', labelsize=20)
        
        plt.savefig(f'eosimg/{base_name}.jpeg', dpi = 1500, bbox_inches='tight')

        with open(f"{base_name}.txt", "w") as f:
            print(f'Unreacted parameters = {popt_unreacted}, omega = {omega_unreacted}', file = f)
            print(f'Reacted parameters = {popt_reacted}, omega = {omega_reacted}', file = f)

        import csv
        with open(f"DATA_unreacted_{base_name}.csv", "w", newline="") as f:
            w = csv.writer(f, delimiter = ',')

            w.writerow(['J_unreacted', 'P_unreacted'])
            w.writerows(zip(Jeval, P_total_unreacted))
        
        with open(f"DATA_reacted_{base_name}.csv", "w", newline="") as f:
            w = csv.writer(f, delimiter = ',')

            w.writerow(['J_reacted', 'P_reacted'])
            w.writerows(zip(Jeval, P_total_reacted))

        #send to cluster
        send = input('Send to cluster? (YES) or (NO): ')
        if (send == 'YES'):
            import subprocess as sub
            sub.run(["scp", f"DATA_unreacted_{base_name}.csv", "gonz1075@bell.rcac.purdue.edu:/scratch/bell/gonz1075/mist/r2d/images"])
            sub.run(["scp", f"DATA_reacted_{base_name}.csv", "gonz1075@bell.rcac.purdue.edu:/scratch/bell/gonz1075/mist/r2d/images"])
        else:
            print('DONE!')
plt.show()

####