import numpy as np
import matplotlib.pyplot as plt



nf = raw_input("Please input range for q,r (nf):")
num_gridpoints = 2001

q = np.linspace(-nf,nf,num_gridpoints)
r = np.linspace(-nf,nf,num_gridpoints)


Roots = np.zeros((num_gridpoints, num_gridpoints))
zeros = 0
for i in xrange(num_gridpoints):
    for j in xrange(num_gridpoints):
        roots = np.roots(np.array([-r[j], q[i]*r[j], -(q[i]+r[j]), q[i]*r[j]]))
        if np.sum(np.abs(np.imag(roots))) > .0001:
            Roots[i,j] = -1.
        elif len(roots) == 1:
            zeros += 1
        else:
            Roots[i,j] = 1.

plt.imshow(Roots.T[:,::-1], cmap='spectral')
plt.colorbar()
plt.savefig('RootsPlot.png')
plt.show()