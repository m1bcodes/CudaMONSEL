import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mim

#data = np.loadtxt("CudaMONSEL\outputs\data54.txt")
#SE = data[:, 4]

data = np.loadtxt("CudaMONSEL\\119-8-13_10-0-56\\output2.txt", max_rows=1)
SE = data

print(SE)
c = 512
r = (int)(SE.shape[0]/c)

plt.plot(SE[1:c])
plt.ylabel('counts')
plt.show()

# SE = SE/np.max(SE)
SE = SE / np.linalg.norm(SE)
SE = SE.reshape((r, c))
plt.imshow(SE, cmap='gray')
plt.show()

mim.imsave("CudaMONSEL\\119-8-13_10-0-56\\output2.png",  SE, cmap='gray')