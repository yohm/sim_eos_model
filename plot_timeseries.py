import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dat = np.loadtxt("timeseries.dat", comments='#', delimiter=' ')

plt.subplot(2,2,1)

plt.xlabel("Time (thousand steps)")
plt.ylabel("Number of species")
plt.plot( dat[:,0]/1000, dat[:,1], 'o-', label='# of species' )

plt.subplot(2,2,2)
plt.xlabel("Time (thousand steps)")
plt.ylabel("<f>")
plt.plot( dat[:,0]/1000, dat[:,2], 'o-', label='<f>' )

plt.subplot(2,2,3)
plt.xlabel("Time (thousand steps)")
plt.ylabel("<k>")
plt.plot( dat[:,0]/1000, dat[:,3], 'o-', label='<k>' )

plt.subplot(2,2,4)
plt.xlabel("Time (thousand steps)")
plt.ylabel("AC")
plt.plot( dat[:,0]/1000, dat[:,4], 'o-', label='AC' )

plt.tight_layout()
plt.savefig("timeseries.png")

