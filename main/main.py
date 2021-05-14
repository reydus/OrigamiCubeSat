import sympy as sp
import sympy.physics.mechanics as mech
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from matplotlib import style
import os.path
from datetime import datetime


class Panel:
    def __init__(self, frame, dimensions=[5], connections=[], state=0):
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

    # If panel has gone beyond the direction of the bus, it's completely in the shadow.
    if get_angle(panel.frame, bus) < - sp.pi /2:
        force_z = 0
    else:
        #force_z = mech.dot(force, bus.y) * sp.cos(get_angle(panel.frame, bus)) ** 2
        #force_z = (sp.cos(get_angle(panel.frame, bus))*force)
        try:
            force_z = (sp.cos(get_angle(panel.frame, bus))*force).args[0][0][1]
        except IndexError:
            force_z = 0 
    # If the light hits on the back of the panel, discount force caused by a state of 1.
    if force_z < 0:
        force_z = force_z / (1+panel.state)  
    #force_z = mech.dot(force, panel.frame.y).evalf()
    torque_zz = - panel.dimensions[0]*3**(1/2)/6 * force_z
    panel.torque = 0*panel.frame.y + 0*panel.frame.x + torque_zz*panel.frame.z
    return force_z

def main():
    theta_zero = 0 * 120/180 * sp.pi#10/180 * sp.pi
    bus = sun.orientnew("bus", "Axis", [angle_to_sun, sun.z])
    order = [bus]

    # MANUAL GENERATION FOR PANELS.
    
    panel_ors["1"] = Panel(
        dimensions=5,
        frame = bus.orientnew("panel"+str(1), "Axis", [theta_zero, bus.z]),
        state = 1
        )
    order.append(panel_ors["1"].frame)
    
    panel_ors["2"] = Panel(
        dimensions=5,
        frame = bus.orientnew("panel"+str(2), "Axis", [theta_zero, bus.z]),
        state = 1
        )
    order.append(panel_ors["2"].frame)
    
    panel_ors["3"] = Panel(
        dimensions=5,
        frame = bus.orientnew("panel"+str(3), "Axis", [theta_zero, bus.z]),
        state = 1
        )
    order.append(panel_ors["3"].frame)
    
    
    upd = 3.6**-1 # hertz
    t= 1/upd
    omega_0 = 0

    
    for panel in panel_ors:
        
        panel = panel_ors[panel]
        c = 0
        alpha = 0 *panel.frame.z
        #theta = sp.acos(mech.dot(panel.frame.x, bus.x).evalf()) * panel.frame.z
        theta = get_angle(panel.frame, bus) * panel.frame.z
        panel.y_axis_theta = []
        panel.y_axis_force = []
        panel.y_axis_force_absolute = []
        panel.y_axis_moment_srp = []
        panel.y_axis_moment_spr = []
        panel.y_axis_moment_dam = []
        panel.y_axis_speed = []
        panel.timeStep = t
        bus_location = 1.5 * sun.x # location in AU referenced on the sun frame.
        omega_0 = 0 * panel.frame.z
        spring_trigger= 0

        panel.hubState = 0
        panel.hubForce = panel_ors["1"].area * calc_srp(bus_location, bus) * (1+panel.hubState) # Simulate hub as a de-activated panel.


        while True:
            f = panel.area * calc_srp(bus_location, bus) * abs(sp.cos(get_angle(panel.frame, bus)))
            if get_angle(panel.frame, bus) < sp.pi/2:
                f = f * (1+panel.state)

            angled_force = calc_torque(panel, f, bus)
            panel.y_axis_force_absolute.append(angled_force)
            
            f = f - panel.hubForce  # Panel moves relative to the hub, so the force must be the difference of the two forces.
            if type(f) != sp.physics.vector.vector.Vector and f == 0:
                f = 0 * bus.y
            angled_force = calc_torque(panel, f, bus)



            spring = get_spring_moment(panel.frame, bus)
            damper = get_joint_damper_monent(omega_0, panel.frame)

            alpha = (panel.torque+spring+damper) / panel.MOI
            d_theta = omega_0 * t + alpha*t*t/2
            theta += d_theta
            c += 1
            omega_0 = omega_0 + alpha*t
            panel.frame.orient(bus, "Axis", [mech.dot(theta, panel.frame.z), bus.z])
            print(str(theta)+", "+str(omega_0)+", "+str(c))
            


            panel.y_axis_theta.append(get_angle(panel.frame, bus))
            panel.y_axis_force.append(angled_force)
            panel.y_axis_moment_spr.append(mech.dot(spring, bus.z).evalf())
            panel.y_axis_moment_srp.append(mech.dot(panel.torque, bus.z).evalf())
            panel.y_axis_moment_dam.append(mech.dot(damper, bus.z).evalf())
            panel.y_axis_speed.append(mech.dot(omega_0, bus.z).evalf())

            if c > 2000:
                print("Iteration limit reached, escaping...")
                panel.iterations = c
                break
            
            if c > 25 and threshold(np.mean(np.abs(np.diff(panel.y_axis_theta[-convergence_mean:]))), 0, thresh=0.0001):
                print("Solution found after "+str(c)+" iterations")
                #print("Force: "+"%.3g" % (2*n*panel.y_axis_force[-1])+"N")
                panel.iterations = c
                break
            '''
            if abs(get_angle(panel.frame, bus)) > sp.pi.evalf()/2:
                print("Exception ocurred: Panel rotated more than 90 degrees.")
                print("Force: "+"%.3g" % (2*n*panel.y_axis_force[-1])+"N")
                break
            '''
    # Print total pressure force felt.
    totalForce = 0
    for i in panel_ors:
        totalForce += panel_ors[i].y_axis_force_absolute[-1]
    
    totalForce += mech.dot(panel_ors["1"].hubForce, bus.y)

    
    print("Total force due to SRP: "+str(totalForce)+" Newtons")


    # Retrieve maximum number of iterations.
    maxIterations = 0
    for i in panel_ors:
        if maxIterations < panel_ors[i].iterations:
            maxIterations = panel_ors[i].iterations

    # Create array with absolute time values from start of simulation to end.
    x_axis = np.linspace(0.0, maxIterations*t, maxIterations)
    
    # Auto-complete for any panel that converges faster than the slowest one.
    for i in panel_ors:
        panel_ors[i].y_axis_theta = auto_complete(panel_ors[i].y_axis_theta, maxIterations)
        panel_ors[i].y_axis_force = auto_complete(panel_ors[i].y_axis_force, maxIterations)
        panel_ors[i].y_axis_moment_srp = auto_complete(panel_ors[i].y_axis_moment_srp, maxIterations)
        panel_ors[i].y_axis_moment_spr = auto_complete(panel_ors[i].y_axis_moment_spr, maxIterations)
        panel_ors[i].y_axis_moment_dam = auto_complete(panel_ors[i].y_axis_moment_dam, maxIterations)
        panel_ors[i].y_axis_speed = auto_complete(panel_ors[i].y_axis_speed, maxIterations)
        panel_ors[i].x_axis = x_axis

    ########## PLOT THETAS
    figtheta, axtheta = plt.subplots()  # Create a figure and an axes.
    for i in panel_ors:
        axtheta.plot(x_axis, panel_ors[i].y_axis_theta, label="Panel no. "+i)
    axtheta.set_xlabel('Time (s)')  # Add an x-label to the axes.
    axtheta.set_ylabel('Panel angle (rads)')  # Add a y-label to the axes.
    axtheta.set_title("Angle plot against time, initial theta ="+str(round(theta_zero.evalf(),2)))  # Add a title to the axes.
    axtheta.grid(True)
    axtheta.legend()

    ########## PLOT FORCES
    figforce, axforce = plt.subplots()

    for i in panel_ors:
        axforce.plot(x_axis, panel_ors[i].y_axis_force, label="Panel no. "+i)
    axforce.set_xlabel("Time (s)")
    axforce.set_ylabel("SRP Force (N)")
    axforce.set_title("Force plot against time, initial theta ="+str(round(theta_zero.evalf(),2)))
    axforce.grid(True)
    axforce.legend()

    ########## PLOT MOMENTS
    figmoment, axmoment = plt.subplots()
    for i in panel_ors:
        axmoment.plot(x_axis, panel_ors[i].y_axis_moment_spr, label="Spring Moment no. "+i)
        axmoment.plot(x_axis, panel_ors[i].y_axis_moment_srp, label="SRP Moment no. "+i)
        axmoment.plot(x_axis, panel_ors[i].y_axis_moment_dam, label="Damper Moment no. "+i)
    axmoment.set_xlabel("Time (s)")
    axmoment.set_ylabel("Moment (N*m)")
    axmoment.set_title("Moment plot against time, initial theta ="+str(round(theta_zero.evalf(),2)))
    axmoment.legend()
    axmoment.grid(True)

    ########## PLOT SPEEDS
    figspeed, axspeed = plt.subplots()
    for i in panel_ors:
        axspeed.plot(x_axis, panel_ors[i].y_axis_speed, label="Angular speed no. "+i)
    axspeed.set_xlabel("Time (s)")
    axspeed.set_ylabel("Angular velocity (Rads/s)")
    axspeed.set_title("Angular velocity plot against time, initial theta ="+str(round(theta_zero.evalf(),2)))
    axspeed.grid(True)

    draw_scene(panel_ors, bus)

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
    c = 0.0277 # for a damp of 0.45, Dampening = C/(2*sqrt(MOI*k)) 
    damper_moment = -c * omega
    return damper_moment

def get_spring_moment(frame1, frame2, alpha=0):
    k = 9.5E-4 #100 * 4.29E-6 # spring constant. N*M/rad
    ref_spring_theta = 0 # rads
    theta = get_angle(frame1, frame2)
    moment = -k * (theta - ref_spring_theta) * frame1.z
    return moment

def get_angle(frame1, frame2):
    angle = sp.acos(mech.dot(frame1.y, frame2.y)).evalf()
    try:
        projection = frame1.y.express(frame2).args[0][0][0] // abs(frame1.y.express(frame2).args[0][0][0]) # If the y axis of the panel has a non-zero projection on the x-axis of the bus, it is at an angle. It falls -ve if CCW, it falls +ve if CW.
    except ZeroDivisionError:
        projection = 1 # if the angle of the panel is perfectly 0, assume CW direction.

    if projection > 0:
        angle = -angle
    return angle

def get_angle_old(frame1, frame2):
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
    import mpl_toolkits.mplot3d.axes3d as p3
    import matplotlib.animation as animation
    import math
    def make_triangles(alpha=90/180 * math.pi, transformation=0):
        #alpha = 90/180 * math.pi
        #print("Painting triangle with "+str(alpha)+"rads")
        tria = np.zeros((3,3))
        h = 3**(1/2)/2 * l

        z_height = -math.sin(alpha) * h  # Height of the piece's ends above XY plane.
        h_xy = math.cos(alpha) * h

        tria[0] = [0, 3**(1/2)/3 * l, 0]
        tria[1] = [(-math.cos(1/6 * math.pi)* h_xy - 1/4 * l),
                    math.sin(1/6 * math.pi) * h_xy + 3**(1/2)/12 * l,
                    z_height]
        tria[2] = [-l/2, - 3**(1/2)/3 * 0.5 * l, 0]

        while transformation != 0:
            theta = 120/180 * math.pi
            c, s = np.cos(theta), np.sin(theta)
            R = np.array(((c, -s, 0), (s, c, 0), (0,0,1)))
            tria[0] = np.dot(R,tria[0])
            tria[1] = np.dot(R,tria[1])
            tria[2] = np.dot(R,tria[2])
            transformation -= 1

        piece = p3.art3d.Poly3DCollection([tria])
        piece.set_color("w")
        piece.set_edgecolor("r")
        return piece

    def update_lines(num, panel_ors):
        while len(ax.collections) > 1:
            ax.collections.pop(1)

        for i in panel_ors:
            piece = make_triangles(panel_ors[i].y_axis_theta[num], int(i)-1)
            if panel_ors[i].state == 0:
                piece.set_color("cornflowerblue")
            elif panel_ors[i].state == 1:
                piece.set_color("blue")
            ax.add_collection3d(piece)
        timeElapsed = "{:.2f}".format(panel_ors["1"].x_axis[num]).zfill(timeStampLen)
        timeLabel.set_text("t = "+timeElapsed+" seconds")
        frameRef = animation_id+"_"+str(str(num).zfill(4))
        fname = "fourth\\frames\\"+frameRef+".png"
        '''
        if not os.path.isfile(fname):
            print("Saving file "+fname)
            fig.savefig(fname)
        '''
        
    # Attaching 3D axis to the figure
    fig = plt.figure()
    ax = p3.Axes3D(fig)

    # hub triangle
    hubcoor = np.zeros((3,3))

    l = panel_ors["1"].dimensions[0]
    hubcoor[0] = [0, 3**(1/2)/3 * l, 0]
    hubcoor[1] = [-l/2, - 3**(1/2)/3 * 0.5 * l, 0]
    hubcoor[2] = [l/2, - 3**(1/2)/3 * 0.5 * l, 0]

    hub = p3.art3d.Poly3DCollection([hubcoor])
    if panel_ors["1"].hubState == 1:
        hub.set_color("blue")
    else:
        hub.set_color("cornflowerblue")
    
    hub.set_edgecolor("k")
    ax.add_collection3d(hub)

    for i in panel_ors:
        piece = make_triangles(panel_ors[i].y_axis_theta[-1], int(i)-1)
        #ax.add_collection3d(piece)
        if panel_ors[i].state == 0:
            piece.set_color("cornflowerblue")
        elif panel_ors[i].state == 1:
            piece.set_color("blue")
        ax.add_collection3d(piece)
    

    # Setting the axes properties
    ax.set_xlim3d([-l, l])
    ax.set_xlabel('X')

    ax.set_ylim3d([-l, l])
    ax.set_ylabel('Y')

    ax.set_zlim3d([-l, l])
    ax.set_zlabel('Z')

    ax.set_title('3D Test')
    timeLabel = ax.text2D(0.05, 0.95, "t = {:.2f}".format(0.00), transform=ax.transAxes)

    
    timeStampLen = len(str(int(panel_ors["1"].x_axis[-1] // 1 + 1))) + 3 # number of digits of last timestamp rounded up to the highest unit, plus one period plus 2 decimal digits.
    #fig.savefig("fourth\\frames\\test frame.png")
    # Creating the Animation object
    line_ani = animation.FuncAnimation(fig, update_lines, panel_ors["1"].iterations, fargs=(panel_ors,),
                                   interval=50, blit=False)

    plt.show()

def auto_complete(array, intendedLength):
    while len(array) < intendedLength:
        array.append(array[-1])
    return array


if __name__ == "__main__":
    #test_get_angle()
    global convergence_mean
    convergence_mean=150
    global sun
    sun = mech.ReferenceFrame("sun")
    global angle_to_sun
    angle_to_sun = 0 #30/180 *sp.pi

    n = 3

    animation_id = [datetime.now().day,datetime.now().month, datetime.now().hour, datetime.now().minute, datetime.now().second]
    
    for i in range(0, len(animation_id)):
        animation_id[i] = str(animation_id[i])
    
    animation_id = "_".join(animation_id)

    
    panel_ors = {}
    panel_ors, bus = main()
    

