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
p_type = data[:, 4]
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

scat_normal = ax.scatter([], [], [], c='blue', marker='o', label='Normal')  
scat_chaotic = ax.scatter([], [], [], c='red', marker='^', label='Chaotic')

def update(k_value):
    global scat_normal, scat_chaotic

    if scat_normal: scat_normal.remove()
    if scat_chaotic: scat_chaotic.remove()


    mask_iter = np.isclose(k, k_value)
    

    mask_normal = mask_iter & (p_type == 0)
    mask_chaotic = mask_iter & (p_type == 1)

    scat_normal = ax.scatter(x[mask_normal], y[mask_normal], z[mask_normal], c='blue', marker='o', s=30, alpha=0.6)

    scat_chaotic = ax.scatter(x[mask_chaotic], y[mask_chaotic], z[mask_chaotic], c='red', marker='^', s=50, alpha=0.8)

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