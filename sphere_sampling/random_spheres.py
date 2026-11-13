# %%
import numpy as np
import matplotlib.pyplot as plt


def cart_rand(max=1):
    coords = np.random.uniform(-max, max, 3)
    coords /= np.linalg.norm(coords)
    return coords


def cart_sphere():
    coords = np.random.normal(0, 1, 3)
    coords /= np.linalg.norm(coords)
    return coords


def polar_rand():
    phi = np.random.uniform(0, 2*np.pi)
    theta = np.random.uniform(0, np.pi)
    coords = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    return coords


def polar_sphere():
    phi = np.random.uniform(0, 2*np.pi)
    theta = np.arccos(np.random.uniform(-1, 1))
    coords = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    return coords


# %%
n = 10000
vectors = np.array([cart_rand() for _ in range(n)])
vectors1 = np.array([cart_sphere() for _ in range(n)])
vectors2 = np.array([polar_sphere() for _ in range(n)])
vectors3 = np.array([polar_rand() for _ in range(n)])

fig = plt.figure(figsize=(40, 20))
ax = fig.add_subplot(411, projection="3d")
ax.scatter(vectors[:, 0], vectors[:, 1], vectors[:, 2], alpha=0.55, s=24)
ax.set(xlabel="x", ylabel="y", zlabel="z", title="Cartesian Random")
ax.set_box_aspect((1, 1, 1))

ax1 = fig.add_subplot(412, projection="3d")
ax1.scatter(vectors1[:, 0], vectors1[:, 1], vectors1[:, 2], alpha=0.55, s=24, c="pink")
ax1.set(xlabel="x", ylabel="y", zlabel="z", title="Sphere-Corrected Cartesian")
ax1.set_box_aspect((1, 1, 1))

ax2 = fig.add_subplot(413, projection="3d")
ax2.scatter(vectors2[:, 0], vectors2[:, 1], vectors2[:, 2], alpha=0.55, s=24, c="pink")
ax2.set(xlabel="x", ylabel="y", zlabel="z", title="Sphere-Correct Polars")
ax2.set_box_aspect((1, 1, 1))

ax3 = fig.add_subplot(414, projection="3d")
ax3.scatter(vectors3[:, 0], vectors3[:, 1], vectors3[:, 2], alpha=0.55, s=24)
ax3.set(xlabel="x", ylabel="y", zlabel="z", title="Sphere-Correct Polars")
ax3.set_box_aspect((1, 1, 1))
plt.tight_layout()
plt.show()
# %%
# Bin by azimuth and z: for unit vectors these bins have equal solid angle.
n_phi, n_z = 24, 12
test_vectors = vectors3
phi = np.mod(np.arctan2(test_vectors[:, 1], test_vectors[:, 0]), 2 * np.pi)
z = np.clip(test_vectors[:, 2], -1, 1)
phi_edges = np.linspace(0, 2 * np.pi, n_phi + 1)
z_edges = np.linspace(-1, 1, n_z + 1)
phi_bin = np.clip(np.digitize(phi, phi_edges) - 1, 0, n_phi - 1)
z_bin = np.clip(np.digitize(z, z_edges) - 1, 0, n_z - 1)
counts, _, _ = np.histogram2d(phi, z, bins=(phi_edges, z_edges))
local_density = counts[phi_bin, z_bin]

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
points = ax.scatter(
    test_vectors[:, 0], test_vectors[:, 1], test_vectors[:, 2],
    c=local_density, cmap="viridis", alpha=0.7, s=24,
)
fig.colorbar(points, ax=ax, label="Points in equal-area angular bin")
ax.set(xlabel="x", ylabel="y", zlabel="z", title="Unit vectors coloured by local density")
ax.set_box_aspect((1, 1, 1))
plt.tight_layout()
plt.show()

# %%
