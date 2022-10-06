import matplotlib
#Prevents the plot from being shown in the screen when saving it to file
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import csv

#Import CSV data
with open("./stats/fractionInfected_gnp_N1000_AG100_T0.03_G0.0006_L1_Wi0.4_Ws0.4_STime50000_R1.csv", "r") as i :
	rawdata = list(csv.reader(i, delimiter = "\t"))

myData = np.array(rawdata[1:], dtype = np.float64)
timeData = myData[:, 1]
infAgSimul = myData[:, 2]
#infSiteSimul = myData[:, 3]
cumSumAg = np.cumsum(infAgSimul)
cumMeanAg = cumSumAg / np.arange(1, len(timeData)+1)

#cumSumSites = np.cumsum(infSiteSimul)
#cumMeanSites = cumSumSites / np.arange(1, len(timeData)+1)

with open("./stats/Runge-Kutta_gnp_N1000_AG100_T0.03_G0.0006_L1_Wi0.4_Ws0.4_STime50000_R1.csv", "r") as j :
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
plt.xlim(0, 50000)
plt.ylim(0, 1)
plt.plot(timeData, infAgSimul, label = "Simulation")
#plt.plot(timeData, infSiteSimul, label = "InfSites")
plt.plot(timeData, cumMeanAg, label = "Cumul. Average")
#plt.plot(timeData, cumMeanSites, label = "Cum.Av.#infSites")
#plt.plot(timeData, infSiteSimul, label = "Model")
#plt.plot(timeData, [np.mean(infAgSimul) for i in range(len(timeData))], label = "Av.#infAg")
plt.plot(timeRK, infAgRK, label = "Model")
#plt.plot(timeRK, infSiteRK, label = "Model-Site")
plt.legend()
plt.grid()
#plt.xscale("log")
#plt.yscale("log")

plt.savefig("./plots/fractionInfected_gnp_N1000_AG100_T0.03_G0.0006_L1_Wi0.4_Ws0.4_STime50000_R1.pdf")
#plt.show()
