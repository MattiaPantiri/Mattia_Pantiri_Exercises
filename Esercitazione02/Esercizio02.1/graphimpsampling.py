import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x, f, error = np.loadtxt("impsampling.txt", usecols=(0,1,2), delimiter='	', unpack='true')
plt.errorbar(x,f,yerr=error)
plt.xlabel('#blocks')
plt.ylabel('I')
plt.grid(True)

plt.show()
