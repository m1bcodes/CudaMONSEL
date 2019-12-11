import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fn1d = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\IIMFPFullPennInterpSiSI.csv"
fn2d = "C:\\Program Files\\NIST\\JMONSEL\\ScatteringTables\\SiTables\\interpNUSimReducedDeltaEFullPennSiSI.csv"

def get_data(fn):
    with open(fn) as f:
        dim = int(f.readline())
        d = np.zeros((dim), dtype=int)
        x = {}
        for i in range(0, dim):
            d[i] = int(f.readline())
            x[i] = np.zeros((d[i]))
            for j in range(0, d[i]):
                x[i][j] = float(f.readline())

        y = np.zeros(d)
        if dim == 1:
            for i in range(0, d[0]):
                y[i] = float(f.readline())
        elif dim == 2:
            for i in range(0, d[0]):
                for j in range(0, d[1]):
                    y[i][j] = float(f.readline())

        return dim, d, x, y

# 1d
dim_, d_, x_, y_ = get_data(fn1d)
plt.plot(x_[0], y_)
plt.show()

# 2d
dim, d, x, y = get_data(fn2d)

x0, x1 = np.meshgrid(x[1], x[0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#ax.scatter(x0, x1, y, cmap='viridis', linewidth=0.5)
ax.contour3D(x0, x1, y, 50, cmap='binary')
#ax.plot_surface(x0, x1, y, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

plt.plot(x[0], np.sum(y, axis=1))
#plt.plot(x_[0], y_)
plt.show()
