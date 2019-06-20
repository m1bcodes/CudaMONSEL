import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mim

data = np.loadtxt("CudaMONSEL\outputs\data25.txt")
SE = data[:, 4]

c = 80
r = (int)(SE.shape[0]/c)

plt.plot(SE[1:c])
plt.ylabel('counts')
plt.show()

# SE = SE/np.max(SE)
SE = SE / np.linalg.norm(SE)
SE = SE.reshape((r, c))
plt.imshow(SE, cmap='gray')
plt.show()

mim.imsave("SE25.png",  SE, cmap='gray')