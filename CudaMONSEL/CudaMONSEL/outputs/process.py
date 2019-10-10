import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mim

#data = np.loadtxt("CudaMONSEL\outputs\data54.txt")
#SE = data[:, 4]

for i in range(0, 100):
    data = np.loadtxt("CudaMONSEL\\outputs\\119-9-9_17-41-50\\output" + str(i) + ".txt", max_rows=1)
    SE = data

    print(SE)
    c = 512
    r = (int)(SE.shape[0]/c)

    # plt.plot(SE[1:c])
    # plt.ylabel('counts')
    # plt.show()

    SE = SE / np.linalg.norm(SE)
    SE = SE.reshape((r, c))
    # plt.imshow(SE, cmap='gray')
    # plt.show()

    mim.imsave("CudaMONSEL\\outputs\\119-9-9_17-41-50\\output" + str(i) + ".png",  SE, cmap='gray')