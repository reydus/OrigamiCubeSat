import sympy as sp
import sympy.physics.mechanics as mech
import matplotlib.pyplot as plt
import numpy as np


bus = mech.ReferenceFrame("bus")


f = 5*bus.y 
n = 1

panel_ors = {}

class Panel:
    def __init__(self, frame, dimensions=[1, 1], connections=[]):
        self.dimensions = dimensions
        self.density = 6 # kg per meter squared
        self.MOI = 1/3 * self.density * dimensions[0]*dimensions[1]
        self.frame = frame
        self.torque = 0



def calc_torque(panel, force):
    force_z = mech.dot(f, panel.frame.y).evalf()
    torque_zz = -panel.dimensions[1]/2 * force_z
    panel.torque = 0*panel.frame.y + 0*panel.frame.x + torque_zz*panel.frame.z
    return force_z
def main():
    for i in range(0, n):
    # for each panel, create a new reference frame.
        panel_ors[i] = Panel(
        dimensions=[1,1],
        connections=[0,"bus",0, 0],
        frame = bus.orientnew(str(i), "Axis", [80/180 * sp.pi, bus.z])
        )
    upd = 100 # hertz
    t= 1/upd
    omega_0 = 0
    
    for panel in panel_ors:
        
        panel = panel_ors[panel]
        iter = 0
        theta = sp.acos(mech.dot(panel.frame.x, bus.x).evalf()) * panel.frame.z
        y_axis_theta = []
        y_axis_force = []
        while mech.dot(theta, bus.z).evalf() > 0:       # When the x-axis of of the bus projects on the y-axis of the bus, the panel latches into place.
            print("current angle is "+str(mech.dot(theta, bus.z).evalf())+"\n")
            angled_force = calc_torque(panel, f)
            alpha = panel.torque/ panel.MOI
            d_theta = omega_0 * t + alpha*t*t/2
            theta += d_theta
            iter += 1
            omega_0 = omega_0 + alpha*t
            panel.frame.orient(bus, "Axis", [mech.dot(theta, panel.frame.z), bus.z])
            print(str(theta)+", "+str(omega_0)+", "+str(iter))
            y_axis_theta.append(theta.magnitude().evalf())
            y_axis_force.append(sp.Abs(angled_force))
        print("Solution found after "+str(iter)+" iterations")
        #x_axis = list(frange(0, t*iter, t))
        x_axis = np.linspace(0.0, iter*t, iter)

        
        
        fig, ax = plt.subplots()  # Create a figure and an axes.
        ax.plot(x_axis, y_axis_theta, label='theta')  # Plot some data on the axes.
        ax.plot(x_axis, y_axis_force, label='force')  # Plot more data on the axes...
        #ax.plot(x_axis, x**3, label='cubic')  # ... and some more.
        ax.set_xlabel('x label')  # Add an x-label to the axes.
        ax.set_ylabel('y label')  # Add a y-label to the axes.
        ax.set_title("Simple Plot")  # Add a title to the axes.
        ax.legend()  # Add a legend.
        plt.show()
            ## rotate panel frame by theta



if __name__ == "__main__":
    main()


