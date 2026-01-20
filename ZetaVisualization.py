import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull

def get_residues():
    return [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

def get_inversion_pairs(k, R):
    pairs = []
    seen = set()
    for z in R:
        o = (k * pow(z, -1, 90)) % 90
        pair = tuple(sorted([z, o]))
        if pair not in seen:
            seen.add(pair)
            pairs.append((z, o))
    return pairs  # 12 pairs for labeling quads

def get_vertices():
    sqrt2 = np.sqrt(2)
    a = (sqrt2 + 4) / 7

    vertices = []

    # Red (axial): (±1, 0, 0) and perms
    for perm in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
        for signs in [1, -1]:
            v = np.array(perm) * signs
            vertices.append(v)
            for cycle in [np.roll(v, 1), np.roll(v, 2)]:
                vertices.append(cycle)

    vertices = []  # Reset, as above is wrong; use the standard set

    # Group 1: permutations of (±sqrt2, 0, 0)
    for i in range(3):
        for sign in [1, -1]:
            v = np.zeros(3)
            v[i] = sign * sqrt2
            vertices.append(v)

    # Group 2: permutations of (±1, ±1, 0)
    for i in range(3):
        for sign1 in [1, -1]:
            for sign2 in [1, -1]:
                v = np.zeros(3)
                v[(i+1)%3] = sign1
                v[(i+2)%3] = sign2
                vertices.append(v)

    # Group 3: all sign combinations of (±a, ±a, ±a)
    for sign1 in [1, -1]:
        for sign2 in [1, -1]:
            for sign3 in [1, -1]:
                v = np.array([sign1, sign2, sign3]) * a
                vertices.append(v)

    vertices = np.unique(vertices, axis=0)  # Remove duplicates if any
    return np.array(vertices)

def visualize_zeta_crystal(k=11, output_file='zeta_crystal_class_11.png'):
    R = get_residues()
    pairs = get_inversion_pairs(k, R)

    vertices = get_vertices()

    # Compute convex hull for surface
    hull = ConvexHull(vertices)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot vertices colored by group
    ax.scatter(vertices[:6,0], vertices[:6,1], vertices[:6,2], color='red', s=50, label='Axial (6)')
    ax.scatter(vertices[6:18,0], vertices[6:18,1], vertices[6:18,2], color='blue', s=50, label='Planar (12)')
    ax.scatter(vertices[18:,0], vertices[18:,1], vertices[18:,2], color='yellow', s=50, label='Octahedral (8)')

    # Plot triangular faces from hull
    for sim in hull.simplices:
        ax.plot_trisurf(vertices[sim,0], vertices[sim,1], vertices[sim,2], alpha=0.2, color='gray')

    # To label faces (triangles, but proxy for quads; label every second for approx 24)
    centroids = []
    for i, sim in enumerate(hull.simplices):
        centroid = np.mean(vertices[sim], axis=0)
        centroids.append(centroid)
        if i < 24:  # Label first 24 triangles with residues
            label = f"{R[i]} ({pairs[i%12][1]})"  # z and o for pair
            ax.text(centroid[0], centroid[1], centroid[2], label, fontsize=8)

    ax.set_title(f'Zeta-Crystal for Class {k} mod 90')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()

    plt.savefig(output_file)
    plt.close(fig)

    print(f"Zeta-Crystal visualization for class {k} saved to {output_file}")

# Example: Generate for class 11
visualize_zeta_crystal(k=11)

# To generate for all 24 classes
for k in get_residues():
    visualize_zeta_crystal(k, output_file=f'zeta_crystal_class_{k}.png')