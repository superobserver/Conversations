import numpy as np
import plotly.graph_objects as go
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
    # Approximate but good enough coordinates for deltoidal icositetrahedron (26 vertices)
    # Based on standard coordinates scaled and permuted
    sqrt2 = np.sqrt(2)
    phi = (1 + np.sqrt(5)) / 2  # golden ratio for aesthetic symmetry

    vertices = []

    # 6 axial vertices: (±c, 0, 0) and permutations
    c = phi
    for perm in [(c, 0, 0), (0, c, 0), (0, 0, c)]:
        for sign in [1, -1]:
            vertices.append(sign * np.array(perm))

    # 12 planar vertices: permutations of (±1, ±1, 0)
    for i in range(3):
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = np.zeros(3)
                v[(i+1)%3] = s1
                v[(i+2)%3] = s2
                vertices.append(v)

    # 8 octahedral-like: (±a, ±a, ±a)
    a = (1 + sqrt2) / 2
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                vertices.append(np.array([s1, s2, s3]) * a)

    return np.array(vertices)

def visualize_zeta_crystal_interactive(k=11, opacity=0.95, output_file='zeta_crystal_class_11.html'):
    R = get_residues()
    pairs = get_inversion_pairs(k, R)

    vertices = get_deltoidal_vertices()
    
    # Compute convex hull for triangulation (will split kites into triangles)
    hull = ConvexHull(vertices)
    
    x, y, z = vertices.T
    i, j, k_tri = hull.simplices.T  # triangulation indices

    # For opacity and gem-like look
    face_colors = []
    for m, tri in enumerate(hull.simplices):
        # Approximate face index (mod 24)
        face_idx = m % 24
        z_val = R[face_idx]
        # Color gradient: darker for small z, brighter for large z
        intensity = (z_val - min(R)) / (max(R) - min(R))
        color = f'rgb({int(100 + 155*intensity)}, {int(120 + 135*intensity)}, {int(180 + 75*intensity)})'
        face_colors.append(color)

    fig = go.Figure()

    # Main opaque crystal mesh
    fig.add_trace(go.Mesh3d(
        x=x, y=y, z=z,
        i=i, j=j, k=k_tri,
        opacity=opacity,           # 0.95 = very opaque, almost solid
        color=face_colors,         # per-face coloring
        colorscale='Viridis',      # fallback if colors don't apply
        flatshading=True,
        lighting=dict(ambient=0.9, diffuse=0.8, specular=0.6, roughness=0.1),
        name='Zeta-Crystal'
    ))

    # Add vertex points for structure
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(size=4, color='darkred', opacity=0.7),
        name='Vertices'
    ))

    # Add hover labels on face centroids (better than fixed text for interactivity)
    centroids = np.mean(vertices[hull.simplices], axis=1)
    hover_labels = []
    for m in range(len(centroids)):
        face_idx = m % 24
        z_val = R[face_idx]
        o_val = (k * pow(z_val, -1, 90)) % 90
        hover_labels.append(f"Face {face_idx+1}<br>z = {z_val}<br>o = {o_val}<br>(class {k})")

    fig.add_trace(go.Scatter3d(
        x=centroids[:,0],
        y=centroids[:,1],
        z=centroids[:,2],
        mode='markers',
        marker=dict(size=0.1, opacity=0),  # invisible points
        hoverinfo='text',
        text=hover_labels,
        name='Face Info'
    ))

    fig.update_layout(
        title=f'Interactive Zeta-Crystal for Class {k} mod 90<br>(Opacity = {opacity} • Rotate & Zoom)',
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='cube',
            bgcolor='black'
        ),
        width=1000,
        height=900,
        showlegend=True,
        paper_bgcolor='black',
        font=dict(color='white')
    )

    # Save as interactive HTML
    fig.write_html(output_file, include_plotlyjs='cdn')
    print(f"Interactive Zeta-Crystal for class {k} saved to {output_file}")
    print("Open the HTML file in any browser → rotate with mouse, zoom with wheel, pan with right-click drag.")

# Generate for class 11 with high opacity (very solid/crystal-like)
visualize_zeta_crystal_interactive(k=11, opacity=0.95)

# To generate for all 24 classes (uncomment if desired):
# for k in get_residues():
#     visualize_zeta_crystal_interactive(k, opacity=0.95, output_file=f'zeta_crystal_class_{k}.html')