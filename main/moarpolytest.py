import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import pylab as pl
import numpy as np




fig = pl.figure()
ax = a3.Axes3D(fig, auto_add_to_figure=False)
ax.set_xlim3d(0, 5)
ax.set_ylim3d(0, 5)
ax.set_zlim3d(0, 5)
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

#fig.add_axes(ax)
#for i in range(3):
vtx = np.random.rand(4,3)
vtx[0] = [1, 1, 2]
vtx[1] = [1, 2, 2]
vtx[2] = [2, 2, 0]
vtx[3] = [2, 1, 0]


tri = a3.art3d.Poly3DCollection([vtx])
tri.set_color(colors.rgb2hex(np.random.rand(3)))
tri.set_edgecolor('k')
ax.add_collection3d(tri)
fig.add_axes(ax)
pl.show()


print("done")
