import numpy as np
import plotly.graph_objects as go
from scipy.spatial import ConvexHull

def get_residues():
    return [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

def get_inversion_for_z(k, z):
    return (k * pow(z, -1, 90)) % 90

def get_deltoidal_vertices():
    """Approximate coordinates for deltoidal icositetrahedron (26 vertices)"""
    sqrt2 = np.sqrt(2)
    phi = (1 + np.sqrt(5)) / 2  # golden ratio

    vertices = []

    # 6 axial: (±phi, 0, 0) and cyclic perms
    c = phi
    for perm in [(c, 0, 0), (0, c, 0), (0, 0, c)]:
        for sign in [1, -1]:
            vertices.append(sign * np.array(perm))

    # 12 planar: permutations of (±1, ±1, 0)
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

    return np.unique(vertices, axis=0)

def visualize_zeta_crystal_interactive(k=11, opacity=0.95, output_file='zeta_crystal_class_11.html'):
    R = get_residues()
    vertices = get_deltoidal_vertices()

    # Convex hull triangulation
    hull = ConvexHull(vertices)
    x, y, z = vertices.T
    i, j, k_tri = hull.simplices.T

    # Per-triangle (face) colors - gradient based on residue value
    face_colors = []
    hover_texts = []
    for m in range(len(hull.simplices)):
        # Cycle through 24 residues
        face_idx = m % 24
        z_val = R[face_idx]
        o_val = get_inversion_for_z(k, z_val)

        # Color: blue-ish gradient, darker for small z, brighter for large
        intensity = (z_val - min(R)) / (max(R) - min(R))
        r = int(60 + 195 * intensity)   # low red → higher
        g = int(100 + 155 * intensity)
        b = int(180 + 75 * intensity)
        color_str = f'rgb({r},{g},{b})'
        face_colors.append(color_str)

        # Hover text
        hover_texts.append(
            f"Face ≈{face_idx+1}<br>"
            f"Residue z = {z_val}<br>"
            f"Inversion o = {o_val}<br>"
            f"Class k = {k} mod 90"
        )

    fig = go.Figure()

    # Main crystal mesh - very opaque
    fig.add_trace(go.Mesh3d(
        x=x, y=y, z=z,
        i=i, j=j, k=k_tri,
        facecolor=face_colors,      # <--- this is the correct way for per-face RGB
        opacity=opacity,            # 0.95 = almost completely solid
        flatshading=True,
        lighting=dict(ambient=0.9, diffuse=0.7, specular=0.6, roughness=0.05),
        name='Zeta-Crystal'
    ))

    # Invisible scatter for hover labels on face centers
    centroids = np.mean(vertices[hull.simplices], axis=1)
    fig.add_trace(go.Scatter3d(
        x=centroids[:,0],
        y=centroids[:,1],
        z=centroids[:,2],
        mode='markers',
        marker=dict(size=1, opacity=0),  # invisible
        hoverinfo='text',
        text=hover_texts,
        name='Face Details (hover)'
    ))

    # Vertex markers
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(size=5, color='darkred', opacity=0.6),
        name='Vertices'
    ))

    fig.update_layout(
        title=f'Zeta-Crystal (Class {k} mod 90) – Interactive 3D<br>(Opacity {opacity} • Rotate / Zoom / Pan)',
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
        font=dict(color='white'),
        hovermode='closest'
    )

    # Save as interactive HTML
    fig.write_html(output_file, include_plotlyjs='cdn')
    print(f"Interactive Zeta-Crystal saved to: {output_file}")
    print("Open in browser → full 3D rotation, zoom, pan. Hover over faces for residue info.")

# Example usage
visualize_zeta_crystal_interactive(k=11, opacity=0.95)

# For all 24 classes (uncomment to generate everything):
# for k in get_residues():
#     visualize_zeta_crystal_interactive(k, opacity=0.95, output_file=f'zeta_crystal_class_{k}.html')