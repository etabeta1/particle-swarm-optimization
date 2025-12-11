import numpy as np
import matplotlib

#matplotlib.use('MacOSX')
#matplotlib.use('TkAgg') for Windows

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = np.loadtxt("points_xy.txt")

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]
k = data[:, 3]
unique_k = np.unique(k)


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

xmin, xmax = np.min(x), np.max(x)
ymin, ymax = np.min(y), np.max(y)
zmin, zmax = np.min(z), np.max(z)


ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin - 1, zmax + 1)  
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")


X_grid = np.linspace(xmin, xmax, 100)
Y_grid = np.linspace(ymin, ymax, 100)
X_surf, Y_surf = np.meshgrid(X_grid, Y_grid)

Z_surf = np.sin(np.sqrt(X_surf ** 2 + Y_surf ** 2))

ax.plot_surface(X_surf, Y_surf, Z_surf, cmap='coolwarm', alpha=0.3)

k0 = unique_k[0]
title_obj = ax.set_title(f"k = {k0}")
scat = ax.scatter([], [], [], c=[])  


def update(k_value):
    global scat

    if scat:
        scat.remove()

    mask = np.isclose(k, k_value)

    scat = ax.scatter(x[mask], y[mask], z[mask], c=z[mask], cmap="viridis", s=40, depthshade=False)

    title_obj.set_text(f"k = {k_value}")

    return []

ani = FuncAnimation(
    fig,
    update,
    frames=unique_k,
    interval=500,
    repeat=True,
    blit=False
)

plt.show()