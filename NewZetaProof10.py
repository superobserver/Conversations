import math
import cmath
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt

# Force R as plain list — this fixes the pow() TypeError
R = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]
R = list(R)  # extra safety


face_colors = [
    f'rgb({int(r*255)}, {int(g*255)}, {int(b*255)})' for r, g, b, _ in plt.cm.tab20(np.linspace(0, 1, 24))
]

def get_operators(k):
    operators = []
    seen_pairs = set()
    for z in R:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in R:
                continue
            pair = tuple(sorted([z, o]))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            z_eff = 91 if z == 1 else z
            o_eff = 91 if o == 1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff * o_eff - k) // 90
            # Always use original z for face color/index
            operators.append((l, m, z_eff, k, z))  # ← face_z = z
            if z != o:
                operators.append((l, m, o_eff, k, z))  # ← still face_z = z (original)
        except ValueError:
            continue
    return operators

def compute_amplitudes_and_traces_for_operator(classes, start_index, total_length, target_z=None):
    segment_start = start_index
    segment_end = start_index + total_length
    
    a = 90
    b = -300
    c = 250 - segment_end
    d = (b**2) - (4 * a * c)
    sol2 = (-b + math.sqrt(d)) / (2 * a) if d >= 0 else (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = int(sol2.real) + 1
    
    list_amp = np.zeros(total_length)
    traces = []  # (x, y0, p, face_z)
    
    all_ops = []
    for k in classes:
        all_ops.extend(get_operators(k))
    
    for x in range(1, new_limit + 1):
        for l_val, m_val, z, primitive, face_z in all_ops:
            if target_z is not None and face_z != target_z:
                continue  # Filter for specific operator z
            y0 = 90 * x * x - l_val * x + m_val
            p = z + 90 * (x - 1)
            if p <= 0 or y0 >= segment_end:
                continue
            if y0 < segment_start:
                diff = segment_start - y0
                n = (diff + p - 1) // p
                current = y0 + n * p
            else:
                current = y0
            traces.append((x, y0, p, face_z))
            while current < segment_end:
                if current >= segment_start:
                    idx = current - segment_start
                    list_amp[idx] += 1
                current += p
    
    return np.array(list_amp), traces

def visualize_operator_cylinder(amplitudes, traces, N_h, title, output_file, merged=False):
    side = math.ceil(math.sqrt(N_h))
    grid = np.zeros((side, side))
    grid.flat[:N_h] = amplitudes / (amplitudes.max() + 1e-8) if amplitudes.max() > 0 else grid
    
    # Cylinder: theta = 2π col / side, z = row, r = 1
    theta = 2 * np.pi * np.arange(side) / side
    z_height = np.arange(side)
    Theta, Z = np.meshgrid(theta, z_height)
    R = np.ones_like(Theta)
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)
    C = grid.T  # Transpose for orientation
    
    fig = go.Figure()
    
    # Amplitude surface
    fig.add_trace(go.Surface(
        x=X, y=Y, z=Z,
        surfacecolor=C,
        colorscale='gray_r',
        showscale=False,
        opacity=0.8,
        name='Amplitude Field'
    ))
    
    # Operator traces as spirals
    for x, y0, p, face_z in traces:
    # Safe lookup: try original R first, then check paired o values
        try:
            face_idx = R.index(face_z)
        except ValueError:
            # If face_z is the paired o, find the corresponding z that inverts to it
            for i, z in enumerate(R):
                o = (k * pow(z, -1, 90)) % 90
                if o == face_z:
                    face_idx = i
                    break
            else:
                face_idx = 0  # fallback to first color if still not found (rare)
    
    
    
    
    
    
    
    
    
    
        #face_idx = np.where(R == face_z)[0][0]
        color = face_colors[face_idx]
        t = np.linspace(0, min(10*p, side*1.5), 200)
        wave_y = y0 + t + 0.08 * side * np.sin(2 * np.pi * t / p)
        row = wave_y // side
        col = wave_y % side
        valid = (row >= 0) & (row < side) & (col >= 0) & (col < side)
        theta_trace = 2 * np.pi * col[valid] / side
        z_trace = row[valid]
        x_trace = np.cos(theta_trace)
        y_trace = np.sin(theta_trace)
        fig.add_trace(go.Scatter3d(
            x=x_trace, y=y_trace, z=z_trace,
            mode='lines',
            line=dict(color=color, width=4),
            name=f'Operator {face_z}',
            opacity=0.7 if merged else 1.0
        ))
    
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z (Height)',
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=1.5)
        ),
        width=800,
        height=800,
        showlegend=merged  # Legend only for merged
    )
    
    fig.write_html(output_file)
    print(f"Interactive cylinder visualization saved to {output_file}")

# Main function to generate individual and merged cylinders
def generate_operator_cylinders(classes, h=20, merged=True):
    N_h = 90 * h**2 - 12 * h + 1
    # Merged
    if merged:
        amplitudes, traces = compute_amplitudes_and_traces_for_operator(classes, 0, N_h)
        visualize_operator_cylinder(amplitudes, traces, N_h, f"Merged 24 Operators (Epoch h={h})", 'merged_cylinder.html', merged=True)
    
    # Individual for each of 24 operators
    for i, z in enumerate(R):
        amplitudes, traces = compute_amplitudes_and_traces_for_operator(classes, 0, N_h, target_z=z)
        visualize_operator_cylinder(amplitudes, traces, N_h, f"Operator {z} Cylinder (Epoch h={h})", f'operator_{z}_cylinder.html', merged=False)

# Run parameters
if __name__ == '__main__':
    classes = [11]  # Example class
    h = 20  # Epoch for visualization
    generate_operator_cylinders(classes, h)