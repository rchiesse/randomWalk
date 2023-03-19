import matplotlib
#Prevents the plot from being shown in the screen when saving it to file
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import csv

with open("./stats/Runge-Kutta_BA-12k-10_N12008_AG1500_T1_G0.002_L1_Wi1_Ws1_STime1000_R1.csv", "r") as j :
	rawRK = list(csv.reader(j, delimiter = "\t"))

rkData = np.array(rawRK[1:], dtype = np.float64)
timeRK = rkData[:, 0]
infAgRK = rkData[:, 1]
#infSiteRK = rkData[:, 2]
#Plot
plt.figure(1, dpi = 120)
plt.title("Fraction of Infected Agents over Time")
plt.xlabel("Time")
plt.ylabel("Infected Fraction")
plt.xlim(0, 1000)
plt.ylim(0, 1)
plt.plot(timeRK, infAgRK, label = "Model")
#plt.plot(timeRK, infSiteRK, label = "Model-Site")
plt.legend()
plt.grid()
#plt.xscale("log")
#plt.yscale("log")

plt.savefig("./plots/fractionInfected_BA-12k-10_N12008_AG1500_T1_G0.002_L1_Wi1_Ws1_STime1000_R1.pdf")
#plt.show()


#AVERAGE DURATION X LAMBDA & AVERAGE DURATION X K
plt.pause(0.001)
plt.clf()
plt.figure(1, dpi = 120)
plt.title("Average Duration over Walk Rate")
plt.xlabel("Walk Rate")
plt.ylabel("Average duration")
plt.xlim(1, 10)
plt.ylim(0, 1000)
plt.xscale("log")
plt.legend()
plt.grid()
plt.savefig("./plots/averages/fractionInfected_BA-12k-10_N12008_AG1500_T1_G0.002_L1_Wi1_Ws1_STime1000_R1.pdf")
