import numpy as np
import matplotlib.pyplot as plt

# --- Make input grid for the surface ---
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x, y)

# Define a function z = f(x, y)
Z = np.sin(np.sqrt(X**2 + Y**2))

# --- Create some points to overlay ---
# Example: random (x, y) in the same range, with z on the surface
x_points = np.random.uniform(-5, 5, 100)
y_points = np.random.uniform(-5, 5, 100)
z_points = np.sin(np.sqrt(x_points**2 + y_points**2))

# --- Plot everything on the SAME axes ---
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Surface
ax.plot_surface(X, Y, Z, alpha=0.7, cmap='viridis')  # alpha to see points better

# Points
ax.scatter(x_points, y_points, z_points, s=20)        # you can tweak s for size

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

plt.show()
