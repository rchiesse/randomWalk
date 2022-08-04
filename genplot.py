import matplotlib
#Prevents the plot from being shown in the screen when saving it to file
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import csv

with open("./stats/Runge-Kutta_gnp_N200_AG40000_T3_G500_L1_STime30_R1.csv", "r") as j :
	rawRK = list(csv.reader(j, delimiter = "\t"))

rkData = np.array(rawRK[1:], dtype = np.float64)
timeRK = rkData[:, 0]
infAgRK = rkData[:, 1]
#infSiteRK = rkData[:, 2]
#Plot
plt.figure(1, dpi = 120)
plt.title("Fraction of Infected Agents/Sites over Time")
plt.xlabel("Time")
plt.ylabel("Infected Fraction")
plt.xlim(0, 30)
plt.ylim(0, 1)
plt.plot(timeRK, infAgRK, label = "Model")
#plt.plot(timeRK, infSiteRK, label = "Model-Site")
plt.legend()
plt.grid()
#plt.xscale("log")
#plt.yscale("log")

plt.savefig("./plots/fractionInfected_gnp_N200_AG40000_T3_G500_L1_STime30_R1.pdf")
#plt.show()
