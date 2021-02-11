import numpy as np
import matplotlib.pyplot as plt

x,y = np.meshgrid(np.linspace(-5,5,10),np.linspace(-5,5,10))

u = y/y * 5
v = x/np.sqrt(x**2 + y**2)
fig, axes = plt.subplots()
#axes.plot(x_axis, y_axis_moment_spr, label="Spring moment")
x,y = np.meshgrid(np.linspace(axes.get_xlim()[0],axes.get_xlim()[1],10),np.linspace(axes.get_ylim()[0],axes.get_ylim()[1],10))

u = y/y * mech.dot()
v = x/np.sqrt(x**2 + y**2)
#fig, axes = plt.subplots()
axes.quiver(x,y,u,v)
plt.show()
