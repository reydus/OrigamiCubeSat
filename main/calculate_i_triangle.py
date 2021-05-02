import math

l = 1

S = l/2
h = S / (3**(1/2)/3)


sect = 3

I = 0

rho = 0.062 # kg per meter squared



for k in range(0, 300):
    sect = k
    for element_no in range(0, sect):
        dx = h/sect
        x = dx*element_no
        dy = S - (3**(1/2)/3)*dx*element_no
        m = 2 * dx * dy * rho
        i = m * (dx*element_no)**2
        #print("Element %d is %fm away from left edge, mass %f, inertia %f" % (element_no, x, m, i))
        I+= i
    print("sections: %d, inertia:%f" % (sect, I))
    I = 0