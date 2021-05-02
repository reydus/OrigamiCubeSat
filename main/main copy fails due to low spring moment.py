import sympy as sp
import sympy.physics.mechanics as mech
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from matplotlib import style

#bus = mech.ReferenceFrame("bus")


class Panel:
    def __init__(self, frame, dimensions=[0.7, 0.7], connections=[], state=0):
        self.state = state
        self.dimensions = [5] # dimensions[0]
        self.area = self.dimensions[0] * 3 * (self.dimensions[0]) / (2 * 3**(1/2)) / 2
        self.density = 0.062 # kg per meter squared
        self.mass = self.density * self.area
        self.MOI = get_MOI_triangle(self.dimensions[0], 0.062)
        self.frame = frame
        self.torque = 0



def get_MOI_triangle(l, rho):
        S = l/2
        h = S / (3**(1/2)/3) # Side divided by tan of 30 deg
        sect = 300
        I = 0
        for element_no in range(0, sect):
            dx = h/sect
            x = dx*element_no
            dy = S - (3**(1/2)/3)*dx*element_no
            m = 2 * dx * dy * rho
            i = m * (dx*element_no)**2
            #print("Element %d is %fm away from left edge, mass %f, inertia %f" % (element_no, x, m, i))
            I += i
        return I


def calc_torque(panel, force, bus):
    force_z = mech.dot(force, bus.y) * sp.cos(get_angle(panel.frame, bus))
    if force_z > 0:
        force_z = force_z * (1+panel.state)  # if the light hits on the front of the panel, use state to determine a coefficient of reflection.
    #force_z = mech.dot(force, panel.frame.y).evalf()
    torque_zz = - panel.dimensions[0]*3**(1/2)/6 * force_z
    panel.torque = 0*panel.frame.y + 0*panel.frame.x + torque_zz*panel.frame.z
    return force_z

def main():
    theta_zero = 30/180 * sp.pi#10/180 * sp.pi
    bus = sun.orientnew("bus", "Axis", [angle_to_sun, sun.z])
    order = [bus]

    # MANUAL GENERATION FOR TWO PANELS.
    '''
    panel_ors["0"] = Panel(
        connections="bus",
        dimensions=[0.7,0.7],
        frame = bus.orientnew("0", "Axis", [theta_zero, bus.z]),
        state = 0
    )
    panel_ors["1"] = Panel(
        connections="0",
        dimensions=[0.7,0.7],
        frame = panel_ors[0].frame.orientnew("1", "Axis", [theta_zero, panel_ors["0"].frame.z]),
        state = 0
    )
    '''
    #PANEL GENERATOR FOR LOOP.

    for i in range(1, n+1):
    # for each panel, create a new reference frame.
        panel_ors[i] = Panel(
        dimensions=[0.7,0.7],
        frame = order[i-1].orientnew("panel"+str(i), "Axis", [theta_zero, bus.z]),
        state = 1
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
        y_axis_moment_dam = []
        y_axis_speed = []
        bus_location = 1.5 * sun.x # location in AU referenced on the sun frame.
        omega_0 = 0 * panel.frame.z
        spring_trigger= 0
        while True: 
            f = panel.area * (1+panel.state)* calc_srp(bus_location, bus)
            angled_force = calc_torque(panel, f, bus)
            #spring = get_spring_moment(panel.frame, bus)
            '''
            if get_angle(panel.frame, bus) < 0:
                if mech.dot(omega_0, panel.frame.z) < 0: # once it reaches a desired velocity ccw, turn off.
                    spring = -panel.torque * 6
                else:
                    spring = 0 *panel.frame.z
            else:
                spring = 0 * panel.frame.z
            '''
            spring = get_spring_moment(panel.frame, bus)
            damper = get_joint_damper_monent(omega_0, panel.frame)

            #c = 4.29E-6 # spring constant. N*M/rad
           
            '''
            ref_spring_theta = 0 # rads
            if get_angle(panel.frame, bus) < 0 and spring_trigger == 0: # start breaking
                time_to_stop = 300 # seconds.
                spring_coeff = - alpha.magnitude()* panel.MOI / time_to_stop
                spring = -spring_coeff * (theta - ref_spring_theta)
                spring_trigger = 1
            elif get_angle(panel.frame, bus) > 0:
                spring_trigger = 0
                spring = 0 * panel.frame.z
            '''

            alpha = (panel.torque+spring+damper) / panel.MOI
            d_theta = omega_0 * t + alpha*t*t/2
            theta += d_theta
            c += 1
            omega_0 = omega_0 + alpha*t
            panel.frame.orient(bus, "Axis", [mech.dot(theta, panel.frame.z), bus.z])
            print(str(theta)+", "+str(omega_0)+", "+str(c))
            


            y_axis_theta.append(get_angle(panel.frame, bus))
            y_axis_force.append(angled_force)
            y_axis_moment_spr.append(mech.dot(spring, bus.z).evalf())
            y_axis_moment_srp.append(mech.dot(panel.torque, bus.z).evalf())
            y_axis_moment_dam.append(mech.dot(damper, bus.z).evalf())
            y_axis_speed.append(mech.dot(omega_0, bus.z).evalf())

            if c > 2000:
                print("Iteration limit reached, escaping...")
                break
            
            if c > 1 and threshold(np.mean(np.abs(np.diff(y_axis_theta[-convergence_mean:]))), 0, thresh=0.0001):
                print("Solution found after "+str(c)+" iterations")
                print("Force: "+"%.3g" % (2*n*y_axis_force[-1])+"N")
                break
            
            if abs(get_angle(panel.frame, bus)) > sp.pi.evalf()/2:
                print("Exception ocurred: Panel rotated more than 90 degrees.")
                print("Force: "+"%.3g" % (2*n*y_axis_force[-1])+"N")
                break
        x_axis = np.linspace(0.0, c*t, c)

        figtheta, axtheta = plt.subplots()  # Create a figure and an axes.
        axtheta.plot(x_axis, y_axis_theta)
        axtheta.plot([0,c*t], [y_axis_theta[-1], y_axis_theta[-1]], "r--")
        axtheta.set_xlabel('Time (s)')  # Add an x-label to the axes.
        axtheta.set_ylabel('Panel angle (rads)')  # Add a y-label to the axes.
        axtheta.set_title("Angle plot against time, initial theta ="+str(round(theta_zero.evalf(),2)))  # Add a title to the axes.
        axtheta.grid(True)
        #figtheta.savefig("theta_plot.fig")
        figtheta.savefig("theta_plot.png")
        #axtheta.grid(axis="both", linestyle='-', linewidth=2)

        figforce, axforce = plt.subplots()
        axforce.plot(x_axis, y_axis_force)
        axforce.plot([0, c*t],[y_axis_force[-1],y_axis_force[-1]], "r--")
        axforce.set_xlabel("Time (s)")
        axforce.set_ylabel("SRP Force (N)")
        axforce.set_title("Force plot against time, initial theta ="+str(round(theta_zero.evalf(),2)))
        axforce.grid(True)
        #axforce.savefig("force_plot.fig")
        figforce.savefig("force_plot.png")
        #axforce.grid(axis="both", linestyle='-', linewidth=2)

        figmoment, axmoment = plt.subplots()
        axmoment.plot(x_axis, y_axis_moment_spr, label="Spring moment")
        axmoment.plot(x_axis, y_axis_moment_srp, label="SRP moment")
        axmoment.plot(x_axis, y_axis_moment_dam, label="Damper moment")
        axmoment.plot([0,c*t],[np.mean(y_axis_moment_srp[-convergence_mean:]),np.mean(y_axis_moment_srp[-convergence_mean:])], "r--")
        axmoment.set_xlabel("Time (s)")
        axmoment.set_ylabel("Moment (N*m)")
        axmoment.set_title("Moment plot against time, initial theta ="+str(round(theta_zero.evalf(),2)))
        axmoment.legend()
        axmoment.grid(True)
        #axmoment.savefig("moment_plot.fig")
        figmoment.savefig("moment_plot.png")
        #axmoment.grid(axis="both", linestyle='-', linewidth=2)
        #
        figspeed, axspeed = plt.subplots()
        axspeed.plot(x_axis, y_axis_speed)
        axspeed.plot([0,c*t],[0,0], "r--")
        axspeed.set_xlabel("Time (s)")
        axspeed.set_ylabel("Angular velocity (Rads/s)")
        axspeed.set_title("Angular velocity plot against time, initial theta ="+str(round(theta_zero.evalf(),2)))
        axspeed.grid(True)
        #axspeed.savefig("speed_plot.fig")
        figspeed.savefig("speed_plot.png")
        #axspeed.grid(axis="both", linestyle='-', linewidth=2)
        draw_scene(panel_ors, bus, y_axis_theta)

        #plt.show()
        return panel_ors, bus


def calc_srp(bus_location, bus):
    alpha = sp.acos(mech.dot(bus_location.normalize(), bus.x)).evalf()
    pressure_ref = 4.56E-6 # n/m^2 at 1 AU
    Rs = bus_location.magnitude()
    pressure = pressure_ref * (1/Rs)**2 * sp.cos(alpha)
    pressure = pressure * bus.y
    return pressure

def threshold(a, b, thresh=0.0001):
    return  a > (b-thresh) and a < (b+thresh)

def get_joint_damper_monent(omega, frame):
    c = 1.57E-4 # for a damp of 0.45, Dampening = C/(2*sqrt(MOI*k))
    damper_moment = -c * omega
    return damper_moment

def get_spring_moment(frame1, frame2, alpha=0):
    k = 4.29E-6 # spring constant. N*M/rad
    ref_spring_theta = 0 # rads
    theta = get_angle(frame1, frame2)
    moment = -k * (theta - ref_spring_theta) * frame1.z
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

def draw_scene(panel_ors, sat, angles):
    '''
    bus = [[[-.215, 0], [.215,0]], [[.215, 0],[.215,-0.93]],[[.215,-0.93],[-.215,-0.93]], [[-.215,-0.93],[-.215, 0]]]
    bus = np.array(bus)

    diagram, axes = plt.subplots()
    for i in bus:
        axes.plot(i.T[0], i.T[1], "k-")

    #maxes = 1.1*np.amax(abs(bus), axis = 0)
    panel_vectors = []
    for index,panel in enumerate(panel_ors):
        panel = panel_ors[panel]
        p1 = [0, sp.cos(get_angle(panel.frame, sat))*panel.dimensions[0]]
        p2 = [0,-sp.sin(get_angle(panel.frame, sat))*panel.dimensions[0]]

        panel_vectors.append([p1, p2])
        if index == 0:
            panel_vectors[index][0][0] += .215
            panel_vectors[index][0][1] += .215
        else:
            panel_vectors[index][0][0] += panel_vectors[index-1][1][0]
            panel_vectors[index][1][0] += panel_vectors[index-1][1][0]
            panel_vectors[index][0][1] += panel_vectors[index-1][1][1]
            panel_vectors[index][1][1] += panel_vectors[index-1][1][1]

        if panel.state == 1:
            color= "b"
        else:
            color= "c"
        
        axes.plot(panel_vectors[index][0],panel_vectors[index][1], color + "-")
        axes.plot(panel_vectors[index][0],panel_vectors[index][1],'ob')
        #plot symmetry
        axes.plot(-1*np.array(panel_vectors[index][0]),panel_vectors[index][1], color + "-")
        axes.plot(-1*np.array(panel_vectors[index][0]),panel_vectors[index][1],"ob")

    axes.axis('equal')  #<-- set the axes to the same scale
    #plt.xlim([-maxes[0],maxes[0]]) #<-- set the x axis limits
    #plt.ylim([-maxes[1],maxes[1]]) #<-- set the y axis limits
    #plt.legend(['V'+str(i+1) for i in range(cols)]) #<-- give a legend
    #plt.grid(b=True, which='major') #<-- plot grid lines

    ###########################

    #axes = plt.gca()

    x,y = np.meshgrid(np.linspace(axes.get_xlim()[0],axes.get_xlim()[1],5),np.linspace(axes.get_ylim()[0],axes.get_ylim()[1],5))
    
    u = y/y * float(sp.sin(angle_to_sun).evalf())
    v = x/x * float(sp.cos(angle_to_sun).evalf())
    #fig, axes = plt.subplots()
    axes.quiver(x,y,u,v)
    axes.set_title("Deflection diagram")
    axes.set_xlabel("X axis (m)")
    axes.set_ylabel("Y axis (m)")
    #diagram.savefig("deflection_diagram.fig")
    diagram.savefig("deflection_diagram.png")
    '''

    import mpl_toolkits.mplot3d.axes3d as p3
    import matplotlib.animation as animation
    import math

    def update_lines(num, angles):
        angle = -angles[num]

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
    line_ani = animation.FuncAnimation(fig, update_lines, len(angles), fargs=(angles,),
                                    interval=50, blit=False)

    plt.show()


if __name__ == "__main__":
    #test_get_angle()
    global convergence_mean
    convergence_mean=300
    global sun
    sun = mech.ReferenceFrame("sun")
    global angle_to_sun
    angle_to_sun = 0 #30/180 *sp.pi

    n = 1

    panel_ors = {}
    panel_ors, bus = main()
    

