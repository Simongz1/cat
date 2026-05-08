import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

degree = float(input('Polynomial Degree for fitting: '))
source = input('Provide absolute path for the image to perform data extraction: ')

#load image
img = mpimg.imread(source)
show = plt.imshow(img)

#extract data from the image using ginput

(minup, minus), = plt.ginput(1)
(maxup, maxus), = plt.ginput(1)
points = plt.ginput(0, 0)
up = np.array([p[0] for p in points])
us = np.array([p[1] for p in points])
plt.close()
#normalize the data
upmax, upmin = 3.5, 0
usmax, usmin = 10, 0

#define a scaling function
def ScaleArray(array, xmin, xmax, ymin, ymax):
    x = np.asarray(array, dtype=float)
    new = (x - xmin) / (xmax - xmin) * (ymax - ymin) + ymin
    return new

#normalize
usnorm = ScaleArray(us, minus, maxus, usmin, usmax)
upnorm = ScaleArray(up, minup, maxup, upmin, upmax)

plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
})

#print normalized data
fig, ax = plt.subplots()
ax.plot(upnorm, usnorm, 'b*', markersize=10, label = 'ReaxFF-lg data')
#ax.grid(True)

#perform linear fitting
fitted = np.polyfit(upnorm, usnorm, degree)
linear = np.polyfit(upnorm, usnorm, 1)
upline = np.linspace(upmin, upmax, 100)

#print the polynomial coefficients
for i in range(len(fitted)):
    print(f'degree: {len(fitted)- i - 1}, coefficient = {fitted[i]}')

ax.plot(upline, np.poly1d(fitted)(upline), 'k-', label = f'{degree:.0f} degree fitting')
ax.plot(upline, np.poly1d(linear)(upline), 'r-', label = '1 degree fitting')
#ax.set_title(f'$u_s$ = {a:.3f} $u_p^3$ + {b:.3f} $u_p^2$ + {c:.3f} $u_p$ + {d:.3f}')
ax.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
ax.set_xlabel('$u_p$ (km/s)', fontsize = 20)
ax.set_ylabel('$u_s$ (km/s)', fontsize = 20)
plt.tight_layout()
plt.savefig('usup_lg.jpeg', dpi = 1500, bbox_inches = 'tight')
plt.show()
