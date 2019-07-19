import matplotlib
import matplotlib.pyplot as plt
import numpy as np

n1, n2, n10, n100 = np.loadtxt("standard.txt", usecols=(0,1,2,3), delimiter='	', unpack='true')

n_bins = 50
n, bins, patches = plt.hist(n1, n_bins, range=(0,1))
plt.xlabel('ciao')
plt.ylabel('prova')
plt.title('Histogram loaded from file!')
plt.grid(True)

plt.figure()
n, bins, patches = plt.hist(n2, n_bins, range=(0,1))
plt.xlabel('ciao')
plt.ylabel('prova')
plt.title('Histogram loaded from file!')
plt.grid(True)

plt.figure()
n, bins, patches = plt.hist(n10, n_bins, range=(0.1, 0.9))
plt.xlabel('ciao')
plt.ylabel('prova')
plt.title('Histogram loaded from file!')
plt.grid(True)

plt.figure()
n, bins, patches = plt.hist(n100, n_bins, range=(0.35, 0.65))
plt.xlabel('ciao')
plt.ylabel('prova')
plt.title('Histogram loaded from file!')
plt.grid(True)


plt.show()
