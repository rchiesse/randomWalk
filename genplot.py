import matplotlib
#Prevents the plot from being shown in the screen when saving it to file
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import csv

with open("./stats/Runge-Kutta_CL_N200_AG15000_T2_G40_L1_STime2_R1.csv", "r") as j :
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
plt.xlim(0, 2)
plt.ylim(0, 1)
plt.plot(timeRK, infAgRK, label = "Model-Ag")
#plt.plot(timeRK, infSiteRK, label = "Model-Site")
plt.legend()
plt.grid()
#plt.xscale("log")
#plt.yscale("log")

plt.savefig("./plots/fractionInfected_CL_N200_AG15000_T2_G40_L1_STime2_R1.pdf")
#plt.show()
