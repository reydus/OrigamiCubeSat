import numpy as np
import matplotlib.pyplot as plt





bus = [[[-.215, 0], [.215,0]], [[.215, 0],[.215,-0.93]],[[.215,-0.93],[-.215,-0.93]], [[-.215,-0.93],[-.215, 0]]]
bus = np.array(bus)

for i in bus:
    plt.plot(i.T[0], i.T[1], "k-")

maxes = 1.1*np.amax(abs(bus), axis = 0)

for index,panel in enumerate(panel_ors):
    for i in range(0, index):
        offset += panel_ors[i].dimensions[0]
    offset += .215    
    p1 = [offset, 0]
    p2 = [offset+panel.dimensions[0], 0]
    

#plt.plot(0,0,'ok') #<-- plot a black point at the origin
plt.axis('equal')  #<-- set the axes to the same scale
#plt.xlim([-maxes[0],maxes[0]]) #<-- set the x axis limits
#plt.ylim([-maxes[1],maxes[1]]) #<-- set the y axis limits
#plt.legend(['V'+str(i+1) for i in range(cols)]) #<-- give a legend
#plt.grid(b=True, which='major') #<-- plot grid lines
plt.show()