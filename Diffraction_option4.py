import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from math import sqrt, ceil

# Operators from your sieve (October 2025 paper)
operators = [
    (72, -1, 17), (72, -1, 91),
    (108, 29, 19), (108, 29, 53),
    (72, 11, 37), (72, 11, 71),
    (18, 0, 73), (18, 0, 89),
    (102, 20, 11), (102, 20, 67),
    (138, 52, 13), (138, 52, 29),
    (102, 28, 31), (102, 28, 47),
    (48, 3, 49), (48, 3, 83),
    (78, 8, 23), (78, 8, 79),
    (132, 45, 7), (132, 45, 41),
    (78, 16, 43), (78, 16, 59),
    (42, 4, 61), (42, 4, 77)
]

def epoch_size(h):
    return 90 * h*h - 12*h + 1

def get_rays(h):
    N = epoch_size(h)
    xmax = int(ceil(sqrt(250 * N / 90))) + 10
    rays = []  # List of (start_y_norm, slope_norm, alpha)
    max_p = 0
    for l, m, z in operators:
        for x in range(1, xmax + 1):
            y = 90 * x*x - l * x + m
            if y >= N: continue
            p = z + 90 * (x - 1)
            if p <= 0: continue
            y_norm = y / N  # Normalize start to [0,1]
            alpha = 1 / (90 * x + 1e-6)  # Fading cancellation power
            slope_norm = p / N  # Normalized slope for vis (rise/run over frame)
            rays.append((y_norm, slope_norm, alpha))
            max_p = max(max_p, p)
    # Re-normalize slopes relative to max_p for better vis
    for i in range(len(rays)):
        y_norm, slope_norm, alpha = rays[i]
        rays[i] = (y_norm, slope_norm * (N / max_p), alpha)
    return rays

def animate_rays(h_min=10, h_max=100, h_step=5, num_frames_per_h=1, ray_length=1.0, fps=5, save=True):
    h_values = list(range(h_min, h_max + 1, h_step))
    all_rays = [get_rays(h) for h in h_values]

    fig, ax = plt.subplots(figsize=(8,8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title(f"Ray Diffraction: h={h_min} (Slopes as periods, opacity ~1/(90x))")
    ax.axis('off')
    lines = []

    def init():
        return lines

    def update(frame):
        h_idx = frame // num_frames_per_h
        subframe = frame % num_frames_per_h  # For smooth intra-h fade if >1
        rays = all_rays[h_idx]
        ax.clear()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(f"Ray Diffraction: h={h_values[h_idx]} (Slopes as periods, opacity ~1/(90x))")
        ax.axis('off')
        for y_norm, slope_norm, alpha in rays:
            # Draw ray: from (0, y_norm) with direction (ray_length, slope_norm * ray_length)
            x_end = ray_length
            y_end = y_norm + slope_norm * ray_length
            # Clip to [0,1] if needed, but allow overflow for visual flow
            ax.plot([0, x_end], [y_norm, y_end], color='black', alpha=alpha * (1 - subframe / num_frames_per_h), lw=0.5)
        return []

    total_frames = len(h_values) * num_frames_per_h
    anim = FuncAnimation(fig, update, frames=total_frames, init_func=init, interval=1000/fps, blit=False)
    if save:
        anim.save(f"ray_evolution_h{h_min}-{h_max}.gif", writer=PillowWriter(fps=fps))
        print(f"Saved animation as ray_evolution_h{h_min}-{h_max}.gif")
    plt.show()

# ==================== RUN IT ====================
if __name__ == "__main__":
    h_min = int(input("Enter min h (e.g. 10): ") or 10)
    h_max = int(input("Enter max h (e.g. 100): ") or 100)
    h_step = int(input("Enter step (e.g. 5): ") or 5)
    num_frames_per_h = int(input("Frames per h (1=instant, 5=smooth fade): ") or 1)
    ray_length = float(input("Ray length (1.0=full frame, >1 for extension): ") or 1.0)
    fps = int(input("FPS (e.g. 5): ") or 5)
    animate_rays(h_min, h_max, h_step, num_frames_per_h, ray_length, fps, save=True)