import plotly.graph_objects as go
import numpy as np

# Residues (24 coprime to 90)
R = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

def get_inversion_pairs(k):
    pairs = []
    seen = set()
    for z in R:
        o = (k * pow(z, -1, 90)) % 90
        pair = tuple(sorted([z, o]))
        if pair not in seen:
            seen.add(pair)
            pairs.append((z, o))
    return pairs

def get_operator_params(k, z):
    # From paper tables for class k=11 (adapt for other k)
    # Example for k=11: (l, m, z) for each z
    # Dummy values; replace with actual from Table 2
    o = (k * pow(z, -1, 90)) % 90
    z_eff = 91 if z == 1 else z
    o_eff = 91 if o == 1 else o
    l = 180 - (z_eff + o_eff)
    m = 90 - (z_eff + o_eff) + (z_eff * o_eff - k) // 90
    return l, m

def visualize_entanglement_configurations(k=11, output_file='zeta_crystal_entanglement.html'):
    pairs = get_inversion_pairs(k)
    
    # 3D positions for nodes (residues) - arrange in a spherical layout for polyhedral feel
    theta = np.linspace(0, 2*np.pi, len(R))
    phi = np.linspace(0, np.pi, len(R))
    r = 1  # Radius
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    
    # Node trace
    node_trace = go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers+text',
        marker=dict(size=8, color='blue', opacity=0.8),
        text=[str(r) for r in R],
        textposition='middle center',
        hoverinfo='text',
        hovertext=[f"Residue: {r}<br>Pair: {get_inversion_pairs(k)[i % len(pairs)]}" for i, r in enumerate(R)]
    )
    
    # Edge traces for entanglement (inversion pairs)
    edge_traces = []
    for z, o in pairs:
        i_z = R.index(z)
        i_o = R.index(o)
        edge_x = [x[i_z], x[i_o]]
        edge_y = [y[i_z], y[i_o]]
        edge_z = [z[i_z], z[i_o]]
        
        # Get operator params for visualization
        l, m = get_operator_params(k, z)
        
        edge_trace = go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(width=6, color='red'),
            hoverinfo='text',
            hovertext=f"Entanglement: {z} ↔ {o}<br>Operator: l={l}, m={m}"
        )
        edge_traces.append(edge_trace)
    
    fig = go.Figure(data=[node_trace] + edge_traces)
    
    fig.update_layout(
        title=f'Zeta-Crystal Entanglement Configurations for Class {k} mod 90',
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='cube'
        ),
        width=1000,
        height=900,
        showlegend=False
    )
    
    fig.write_html(output_file)
    print(f"Interactive visualization saved to {output_file}. Open in browser to explore.")

# Generate for class 11
visualize_entanglement_configurations() 