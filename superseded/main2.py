import sympy as sp
import sympy.physics.mechanics as mech
import matplotlib.pyplot as plt
import numpy as np


#bus = mech.ReferenceFrame("bus")
sun = mech.ReferenceFrame("sun")
bus = sun.orientnew("bus", "Axis", [0.25*sp.pi, sun.z])


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
        self.area = dimensions[0]* dimensions[1]




def calc_torque(panel, force):
    force_z = mech.dot(force, panel.frame.y).evalf()
    torque_zz = -panel.dimensions[1]/2 * force_z
    panel.torque = 0*panel.frame.y + 0*panel.frame.x + torque_zz*panel.frame.z
    return force_z

def main():
    bus.orient(sun, "Axis", [30/180 *sp.pi, sun.z])
    for i in range(0, n):
    # for each panel, create a new reference frame.
        panel_ors[i] = Panel(
        dimensions=[1,1],
        connections=[0,"bus",0, 0],
        frame = bus.orientnew("panel"+str(i), "Axis", [10/180 * sp.pi, bus.z])
        )
    upd = 3.6**-1 # hertz
    t= 1/upd
    omega_0 = 0
    

    for panel in panel_ors:
        
        panel = panel_ors[panel]
        c = 0
        theta = sp.acos(mech.dot(panel.frame.x, bus.x).evalf()) * panel.frame.z
        y_axis_theta = []
        y_axis_force = []
        y_axis_moment_srp = []
        y_axis_moment_spr = []
        bus_location = 1.5 * sun.x # location in AU referenced on the sun frame.
        
        
        while True: # get_angle(panel.frame, bus) > -0.18:       # When the x-axis of of the bus projects on the y-axis of the bus, the panel latches into place.
            print("current angle is "+str(mech.dot(theta, bus.z).evalf())+"\n")
            f = panel.area * calc_srp(bus_location, bus, state=0)
            angled_force = calc_torque(panel, f)
            spring = get_spring_moment(panel.frame, bus)
            alpha = (panel.torque+spring) / panel.MOI
            d_theta = omega_0 * t + alpha*t*t/2
            theta += d_theta
            c += 1
            omega_0 = omega_0 + alpha*t
            panel.frame.orient(bus, "Axis", [mech.dot(theta, panel.frame.z), bus.z])
            print(str(theta)+", "+str(omega_0)+", "+str(c))
            y_axis_theta.append(theta.magnitude().evalf())
            y_axis_force.append(angled_force)
            y_axis_moment_spr.append(spring.magnitude().evalf())
            y_axis_moment_srp.append(panel.torque.magnitude().evalf())
        print("Solution found after "+str(c)+" iterations")
        #x_axis = list(frange(0, t*c, t))
        x_axis = np.linspace(0.0, c*t, c)

        
        
        fig, ax = plt.subplots()  # Create a figure and an axes.
        ax.plot(x_axis, y_axis_theta, label='theta')  # Plot some data on the axes.
        ax.plot(x_axis, y_axis_force, label='force')  # Plot more data on the axes...
        ax.plot(x_axis, y_axis_moment_spr, label='moment_spr')
        ax.plot(x_axis, y_axis_moment_srp, label='moment_srp')
        #ax.plot(x_axis, x**3, label='cubic')  # ... and some more.
        ax.set_xlabel('x label')  # Add an x-label to the axes.
        ax.set_ylabel('y label')  # Add a y-label to the axes.
        ax.set_title("Simple Plot")  # Add a title to the axes.
        ax.legend()  # Add a legend.
        plt.show()
            ## rotate panel frame by theta
#def time_to_settle(points, theta, ):
'''
def plotter(x, y):
    for i in y
        plt.plot(x,i)
    plt.show()
'''
def calc_srp(bus_location, bus, state=0):
    alpha = sp.acos(mech.dot(bus_location.normalize(), bus.x)).evalf()
    pressure_ref = 4.56E-6 # n/m^2 at 1 AU
    Rs = bus_location.magnitude()
    pressure = pressure_ref * (1+state) * (1/Rs)**2 * sp.cos(alpha)**2
    pressure = pressure * bus.y
    return pressure

def threshold(a, b, thresh=0.0001):
    return  a > (b-thresh) and a < (b+thresh)

def get_spring_moment(frame1, frame2):
    c = 4.29E-6 # spring constant. N*M/rad
    ref_spring_theta = 0 # rads
    theta = get_angle(frame1, frame2)
    moment = -c * (theta - ref_spring_theta) * frame1.z
    return moment


def get_angle(frame1, frame2):
    angle = sp.asin(mech.dot(frame1.y.express(frame2), frame2.x).evalf())
    sector = mech.dot(frame1.x, frame2.x).evalf()
    sign = sector / abs(sector)
    if sign == -1.0:
        angle = sp.pi - angle
    angle = -angle.evalf()
    if threshold(angle, -sp.pi):
        angle = sp.pi
    return angle # positive when the y axis turns anti-CW from frame2's y-axis.

def test_get_angle():
    bus = mech.ReferenceFrame("bus")
    one = [bus.orientnew("one", "Axis", [-30/180 *sp.pi, bus.z]), -30/180 *sp.pi]
    three = [bus.orientnew("three", "Axis", [-90/180 *sp.pi, bus.z]), -90/180 *sp.pi]
    five = [bus.orientnew("five", "Axis", [-120/180 * sp.pi, bus.z]), -120/180 * sp.pi]
    six_cw = [bus.orientnew("six_cw", "Axis", [-180/180 * sp.pi, bus.z]), 180/180 *sp.pi]
    six_ccw = [bus.orientnew("six_ccw", "Axis", [180/180 * sp.pi, bus.z]), 180/180 *sp.pi]
    seven = [bus.orientnew("seven", "Axis", [150/180 *sp.pi, bus.z]), 150/180 * sp.pi]
    nine = [bus.orientnew("nine", "Axis", [90/180 * sp.pi, bus.z]), 90/180 *sp.pi]
    eleven = [bus.orientnew("eleven", "Axis", [30/180 *sp.pi, bus.z]), 30/180 *sp.pi]

    tests = [one, three, six_cw, six_ccw, nine, eleven]
    for i in tests:
        angle = get_angle(i[0], bus)
        if threshold(angle, i[1]):
            pass
        else:
            print("Angle test did not pass: \n")
            print(i[0].name)
            print("Got: "+str(angle)+", expected: "+str(i[1]))
            exit()
    print("All angle tests for the reference system passed.")
if __name__ == "__main__":
    test_get_angle()

    main()


