import matplotlib
#Prevents the plot from being shown in the screen when saving it to file
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import csv

plt.figure(1, dpi = 120)
plt.title("Fraction of Infected Agents over Time")
plt.xlabel("Time")
plt.ylabel("Infected Fraction")
plt.xlim(0, 5000)
plt.ylim(0, 1)

def incluirPlot(nomeArq, ini, lbl, style) :
	with open(nomeArq, "r") as i :
 		rawdata = list(csv.reader(i, delimiter = "	"))

	myData = np.array(rawdata[1:], dtype = np.float64)
	xData = myData[:, ini]
	yData = myData[:, ini + 1]

	plt.plot(xData, yData, label = lbl, linewidth = (1 + (1 - ini)*0.5), linestyle = style) 


#Import CSV data
#incluirPlot("./stats/fractionInfected_GNP-100k-20_N100000_AG5000_T1_G0.005_L1_Wi0.5_Ws0.5_STime5000_R1.csv", 1, "w0.5", "solid")
incluirPlot("./stats/Runge-Kutta_BLOCK_GNP-100k-20_N100000_AG5000_T1_G0.005_L1_Wi0.5_Ws0.5_STime5000_R1.csv", 0, "w0.5 block", "solid")
incluirPlot("./stats/Runge-Kutta_MASTER_GNP-100k-20_N100000_AG5000_T1_G0.005_L1_Wi0.5_Ws0.5_STime5000_R1.csv", 0, "w0.5 master", "dashed")


plt.legend(fontsize="7")
#plt.legend(fontsize="7", ncol=4, loc = (0.15, 0.15))
plt.grid()
#plt.xscale("log")
#plt.yscale("log")

plt.savefig("./plots/fractionInfected_GNP-100k-20_N100000_AG5000_T1_G0.005_L1_Wi0.5_Ws0.5_STime5000_R1.pdf")
#plt.show()


#AVERAGE DURATION X LAMBDA & AVERAGE DURATION X K
plt.pause(0.001)
plt.clf()
plt.figure(1, dpi = 120)
plt.title("Average Duration over Walk Rate")
plt.xlabel("Walk Rate")
plt.ylabel("Average duration")
plt.xlim(1, 10)
plt.ylim(0, 5000)
plt.xscale("log")
plt.legend()
plt.grid()
plt.savefig("./plots/averages/fractionInfected_GNP-100k-20_N100000_AG5000_T1_G0.005_L1_Wi0.5_Ws0.5_STime5000_R1.pdf")
