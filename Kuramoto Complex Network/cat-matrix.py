import numpy as np
import matplotlib.pyplot as plt

M = np.loadtxt("cat-cortex.csv",unpack=True)

plt.matshow(M,cmap=plt.cm.YlOrRd);
plt.colorbar()
plt.show()

