
def get_force(t):
    if t > 1:
        return 1.52E-6
    else:
        return 0


theta = 0
past = 

while True:
    moment = cos(theta) * get_force(t) * 0.5
    spring = 4.29E-6 * theta
    theta = 


import numpy as np
import matplotlib.pyplot as plt

plt.axis([0, 100, 0, 1])

for i in range(100):
    y = np.random.random()
    plt.scatter(i, y)
    plt.pause(0.05)
    plt.

plt.show()