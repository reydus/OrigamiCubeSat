import numpy as np
import matplotlib.pyplot as plt


gsc = 1362 # kW/m^2

p_srp = gsc / 297000000 # divide by c

planets_radius = [0.39, 0.72, 1, 1.52, 5.2, 9.54, 19.18, 30.06]

planets_names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]

planets_colors = ["darkgray", "orange", "blue", "red", "brown", "moccasin", "teal", "royalblue"]

x = np.linspace(0.2, planets_radius[-1],1000)

y = p_srp * 1 / (x**2)

plt.plot(x, y)

for i in range(0, len(planets_names[:5])):
    plt.axvline(x=planets_radius[i], label=planets_names[i]+" orbit", c=planets_colors[i])


plt.yscale("log")
plt.grid()
plt.ylabel("SRP Pressure (N/m^2)")
plt.xlabel("Orbit radius (AU)")
plt.xlim([0, planets_radius[4]*1.15])
plt.legend()
plt.title("Decrease of SRP as function of orbit radius")
plt.show()