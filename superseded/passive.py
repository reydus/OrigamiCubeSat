import sympy as sp
import sympy.physics.mechanics as mech
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from matplotlib import style



class Panel:
    def __init__(self, frame, dimensions=[0.7, 0.7], connections=[]):
        self.dimensions = dimensions
        self.density = 0.062 # kg per meter squared
        self.MOI = 1/3 * self.density * dimensions[0]*dimensions[1]
        self.frame = frame
        self.torque = 0
        self.area = dimensions[0]* dimensions[1]




def calc_torque(panel, force):
    #force_z = mech.dot(force, bus.y)
    force_z = mech.dot(force, panel.frame.y).evalf()
    torque_zz = -panel.dimensions[1]/2 * force_z
    panel.torque = 0*panel.frame.y + 0*panel.frame.x + torque_zz*panel.frame.z
    return force_z

def main():
    theta_zero = 10/180 * sp.pi
    bus = sun.orientnew("bus", "Axis", [angle_to_sun, sun.z])
    order = [bus]

    for i in range(1, n+1):
    # for each panel, create a new reference frame.
        panel_ors[i] = Panel(
        dimensions=[0.7,0.7],
        frame = order[i-1].orientnew("panel"+str(i), "Axis", [theta_zero, bus.z])
        )
        order.append(panel_ors[i].frame)


    upd = 3.6**-1 # hertz
    t= 1/upd
    omega_0 = 0
    #return panel_ors
    
    for panel in panel_ors:
        
        panel = panel_ors[panel]
        c = 0
        alpha = 0 *panel.frame.z
        theta = sp.acos(mech.dot(panel.frame.x, bus.x).evalf()) * panel.frame.z
        y_axis_theta = []
        y_axis_force = []
        y_axis_moment_srp = []
        y_axis_moment_spr = []
        y_axis_speed = []
        bus_location = 1.5 * sun.x # location in AU referenced on the sun frame.
        
        spring_trigger= 0

        while True: 
            f = panel.area * calc_srp(bus_location, bus, state=0)
            angled_force = calc_torque(panel, f)
            #spring = get_spring_moment(panel.frame, bus)
            spring = -panel.torque * 6
            #c = 4.29E-6 # spring constant. N*M/rad
           


            alpha = (panel.torque+spring) / panel.MOI
            d_theta = omega_0 * t + alpha*t*t/2
            theta += d_theta
            c += 1
            omega_0 = omega_0 + alpha*t
            panel.frame.orient(bus, "Axis", [mech.dot(theta, panel.frame.z), bus.z])
            print(str(theta)+", "+str(omega_0)+", "+str(c))
            


            y_axis_theta.append(get_angle(panel.frame, bus))
            y_axis_force.append(angled_force)
            y_axis_moment_spr.append(spring.magnitude().evalf())
            y_axis_moment_srp.append(panel.torque.magnitude().evalf())
            y_axis_speed.append(mech.dot(omega_0, bus.z).evalf())

            if c > 2000:
                print("Iteration limit reached, escaping...")
                break
            if threshold(np.mean(y_axis_theta[-convergence_mean:-1]), 0, thresh=0.001):
                print("Solution found after "+str(c)+" iterations")
                break
        
        x_axis = np.linspace(0.0, c*t, c)

        figtheta, axtheta = plt.subplots()  # Create a figure and an axes.
        axtheta.plot(x_axis, y_axis_theta)
        axtheta.plot([0,c*t], [0,0], "r--")
        axtheta.set_xlabel('Time (s)')  # Add an x-label to the axes.
        axtheta.set_ylabel('Panel angle (rads)')  # Add a y-label to the axes.
        axtheta.set_title("Angle plot against time, initial theta ="+str(theta_zero))  # Add a title to the axes.

        figforce, axforce = plt.subplots()
        axforce.plot(x_axis, y_axis_force)
        axforce.plot([0, c*t],[y_axis_force[-1],y_axis_force[-1]], "r--")
        axforce.set_xlabel("Time (s)")
        axforce.set_ylabel("SRP Force (N)")
        axforce.set_title("Force plot against time, initial theta ="+str(theta_zero))

        figmoment, axmoment = plt.subplots()
        axmoment.plot(x_axis, y_axis_moment_spr, label="Spring moment")
        axmoment.plot(x_axis, y_axis_moment_srp, label="SRP moment")
        axmoment.plot([0,c*t],[np.mean(y_axis_moment_srp[-convergence_mean:-1]),np.mean(y_axis_moment_srp[-convergence_mean:-1])], "r--")
        axmoment.set_xlabel("Time (s)")
        axmoment.set_ylabel("Moment (N*m)")
        axmoment.set_title("Moment plot against time, initial theta ="+str(theta_zero))
        axmoment.legend()

        figspeed, axspeed = plt.subplots()
        axspeed.plot(x_axis, y_axis_speed)
        axspeed.plot([0,c*t],[0,0], "r--")
        axspeed.set_xlabel("Time (s)")
        axspeed.set_ylabel("Angular velocity (Rads/s)")
        axspeed.set_title("Angular velocity plot against time, initial theta ="+str(theta_zero))

        plt.show()
        return panel_ors, bus


def calc_srp(bus_location, bus, state=0):
    alpha = sp.acos(mech.dot(bus_location.normalize(), bus.x)).evalf()
    pressure_ref = 4.56E-6 # n/m^2 at 1 AU
    Rs = bus_location.magnitude()
    pressure = pressure_ref * (1+state) * (1/Rs)**2 * sp.cos(alpha)**2
    pressure = pressure * bus.y
    return pressure

def threshold(a, b, thresh=0.0001):
    return  a > (b-thresh) and a < (b+thresh)

def get_spring_moment(frame1, frame2, alpha=0):
    c = 4.29E-6 # spring constant. N*M/rad
    ref_spring_theta = 0 # rads
    theta = get_angle(frame1, frame2)
    if theta < 0:
        moment = -c * (theta - ref_spring_theta) * frame1.z
    else:
        moment = 0 * frame1.z
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
'''
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
'''
def draw_scene(panel_ors, sat):
    bus = [[[-.215, 0], [.215,0]], [[.215, 0],[.215,-0.93]],[[.215,-0.93],[-.215,-0.93]], [[-.215,-0.93],[-.215, 0]]]
    bus = np.array(bus)

    for i in bus:
        plt.plot(i.T[0], i.T[1], "k-")

    #maxes = 1.1*np.amax(abs(bus), axis = 0)
    panel_vectors = []
    for index,panel in enumerate(panel_ors):
        panel = panel_ors[panel]
        p1 = [0, sp.cos(get_angle(panel.frame, sat))*panel.dimensions[0]]
        #p2 = [sp.cos(get_angle(panel.frame, sat))*panel.dimensions[0],sp.sin(get_angle(panel.frame, sat))*panel.dimensions[0]]
        p2 = [0,-sp.sin(get_angle(panel.frame, sat))*panel.dimensions[0]]

        panel_vectors.append([p1, p2])
    for i in range(0, len(panel_vectors)):
        if i == 0:
            panel_vectors[i][0][0] += .215
            panel_vectors[i][0][1] += .215
        else:
            panel_vectors[i][0][0] += panel_vectors[i-1][1][0]
            panel_vectors[i][1][0] += panel_vectors[i-1][1][0]
            panel_vectors[i][0][1] += panel_vectors[i-1][1][1]
            panel_vectors[i][1][1] += panel_vectors[i-1][1][1]


    for i in panel_vectors:
        plt.plot(i[0], i[1], "b-")
        plt.plot(i[0],i[1],'ob')
        #plot symmetry
        plt.plot(-1*np.array(i[0]),np.array(i[1]), "b-")
        plt.plot(-1*np.array(i[0]),np.array(i[1]),"ob")
    plt.axis('equal')  #<-- set the axes to the same scale
    #plt.xlim([-maxes[0],maxes[0]]) #<-- set the x axis limits
    #plt.ylim([-maxes[1],maxes[1]]) #<-- set the y axis limits
    #plt.legend(['V'+str(i+1) for i in range(cols)]) #<-- give a legend
    #plt.grid(b=True, which='major') #<-- plot grid lines

    ###########################

    axes = plt.gca()

    x,y = np.meshgrid(np.linspace(axes.get_xlim()[0],axes.get_xlim()[1],5),np.linspace(axes.get_ylim()[0],axes.get_ylim()[1],5))

    u = y/y * float(sp.sin(angle_to_sun).evalf())
    v = x/x * float(sp.cos(angle_to_sun).evalf())
    #fig, axes = plt.subplots()
    axes.quiver(x,y,u,v)
    plt.show()

if __name__ == "__main__":
    #test_get_angle()
    global convergence_mean
    convergence_mean=300
    global sun
    sun = mech.ReferenceFrame("sun")
    global angle_to_sun
    angle_to_sun = 30/180 *sp.pi

    n = 1

    panel_ors = {}
    panel_ors, bus = main()
    draw_scene(panel_ors, bus)

