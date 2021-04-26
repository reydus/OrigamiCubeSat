import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import math

def update_lines(num):
    angle = -1/2 * math.pi + (1/25 * math.pi * num)
    while len(ax.collections) > 1:
        ax.collections.pop(1)
    pieces = make_triangles(angle)
    for i in pieces:
        ax.add_collection3d(i)
    #return lines

# Attaching 3D axis to the figureZ
fig = plt.figure()
ax = p3.Axes3D(fig)

# hub triangle
hubcoor = np.zeros((3,3))

l = 1
hubcoor[0] = [0, 3**(1/2)/3 * l, 0]
hubcoor[1] = [-l/2, - 3**(1/2)/3 * 0.5 * l, 0]
hubcoor[2] = [l/2, - 3**(1/2)/3 * 0.5 * l, 0]

hub = p3.art3d.Poly3DCollection([hubcoor])
hub.set_color("w")
hub.set_edgecolor("k")

def make_triangles(alpha=90/180 * math.pi):
    # first piece
    tria = np.zeros((3,3))
    #alpha =  90/180 * math.pi
    h = 3**(1/2)/2 * l

    z_height = math.sin(alpha) * h  # Height of the piece's ends above XY plane.
    h_xy = math.cos(alpha) * h

    tria[0] = [0, 3**(1/2)/3 * l, 0]
    tria[1] = [-math.cos(1/6 * math.pi)* h_xy + -1/4 * l,
                math.sin(1/6 * math.pi) * h_xy + 3**(1/2)/12 * l,
                z_height]
    tria[2] = [-l/2, - 3**(1/2)/3 * 0.5 * l, 0]

    piece = p3.art3d.Poly3DCollection([tria])
    piece.set_color("w")
    piece.set_edgecolor("r")

    # transformation matrix for 120 degree counter-clockwise replication
    theta = 120/180 * math.pi
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s, 0), (s, c, 0), (0,0,1)))

    # triangle 2 and 3
    tria2 = np.random.rand(3,3)
    tria2[0] = np.dot(R,tria[0])
    tria2[1] = np.dot(R,tria[1])
    tria2[2] = np.dot(R,tria[2])

    piece2 = p3.art3d.Poly3DCollection([tria2])
    piece2.set_color("w")
    piece2.set_edgecolor("r")

    tria3 = np.random.rand(3,3)
    tria3[0] = np.dot(R,tria2[0])
    tria3[1] = np.dot(R,tria2[1])
    tria3[2] = np.dot(R,tria2[2])

    piece3 = p3.art3d.Poly3DCollection([tria3])
    piece3.set_color("w")
    piece3.set_edgecolor("r")

    return [piece, piece2, piece3]

hubs = ax.add_collection3d(hub)

# Setting the axes properties
ax.set_xlim3d([-1, 1])
ax.set_xlabel('X')

ax.set_ylim3d([-1, 1])
ax.set_ylabel('Y')

ax.set_zlim3d([-1, 1])
ax.set_zlabel('Z')

ax.set_title('3D Test')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, update_lines, 25,
                                   interval=50, blit=False)

plt.show()