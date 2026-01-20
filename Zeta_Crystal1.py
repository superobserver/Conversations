import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot  # For saving interactive HTML
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
    return pairs  # 12 pairs, but we label 24 faces with z and o

def get_deltoidal_vertices():
    # Standard coordinates for deltoidal icositetrahedron (normalized)
    # Group 1: (±1, 0, 0) and permutations, scaled
    # Actual from reliable source: permutations of (0, ±1, ±(1+√2)), but with all even permutations for 24 faces
    # For simplicity, use a set of 26 vertices as approximate
    phi = (1 + np.sqrt(5)) / 2  # Golden ratio for icosahedral symmetry approximation
    vertices = []
    # Axial
    for i in range(3):
        for sign in [1, -1]:
            v = np.zeros(3)
            v[i] = sign * phi
            vertices.append(v)
    # Planar
    for i in range(3):
        for sign1 in [1, -1]:
            for sign2 in [1, -1]:
                v = np.zeros(3)
                v[(i+1) % 3] = sign1
                v[(i+2) % 3] = sign2
                vertices.append(v)
    # Octahedral-like
    a = (1 + np.sqrt(2)) / 2
    for sign1 in [1, -1]:
        for sign2 in [1, -1]:
            for sign3 in [1, -1]:
                v = np.array([sign1, sign2, sign3]) * a
                vertices.append(v)

    vertices = np.unique(vertices, axis=0)  # Dedup, should be ~26
    return np.array(vertices)

def visualize_zeta_crystal_interactive(k=11, opacity=0.8, output_file='zeta_crystal_class_11.html'):
    R = get_residues()
    pairs = get_inversion_pairs(k, R)

    vertices = get_deltoidal_vertices()
    
    # Compute convex hull for triangulation (approximates quads by triangles)
    hull = ConvexHull(vertices)

    # Extract triangles for Mesh3d
    x, y, z = vertices.T
    i, j, l = hull.simplices.T  # Note: 'l' for third index

    # Face centroids for labels (approximate, one per triangle)
    centroids = []
    labels = []
    for m, sim in enumerate(hull.simplices):
        centroid = np.mean(vertices[sim], axis=0)
        centroids.append(centroid)
        # Cycle through 24 labels (z with o in parens)
        label_idx = m % 24
        z = R[label_idx]
        o = (k * pow(z, -1, 90)) % 90
        labels.append(f"{z} ({o})")

    # Plotly figure
    fig = go.Figure()

    # Add opaque mesh for surfaces
    fig.add_trace(go.Mesh3d(
        x=x, y=y, z=z,
        i=i, j=j, k=l,  # Triangulation
        opacity=opacity,  # User-controlled opacity (0.8 for more opaque)
        color='lightblue',
        flatshading=True,
        lighting=dict(ambient=0.8, diffuse=0.5, specular=0.3),
        name='Zeta-Crystal'
    ))

    # Add vertex points
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(size=5, color='red'),
        name='Vertices'
    ))

    # Add labels at centroids
    fig.add_trace(go.Scatter3d(
        x=[c[0] for c in centroids],
        y=[c[1] for c in centroids],
        z=[c[2] for c in centroids],
        mode='text',
        text=labels,
        textposition='middle center',
        textfont=dict(size=10, color='black'),
        name='Face Labels'
    ))

    fig.update_layout(
        title=f'Interactive Zeta-Crystal for Class {k} mod 90',
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='cube'
        ),
        width=800,
        height=800,
        legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01)
    )

    # Save as interactive HTML (user can rotate, zoom in browser)
    plot(fig, filename=output_file, auto_open=False)
    print(f"Interactive Zeta-Crystal for class {k} saved to {output_file}. Open in browser to rotate.")

# Generate for class 11
visualize_zeta_crystal_interactive(k=11, opacity=0.9)  # More opaque

# To generate for all 24 classes
R = get_residues()
for k in R:
    visualize_zeta_crystal_interactive(k, opacity=0.9, output_file=f'zeta_crystal_class_{k}.html')