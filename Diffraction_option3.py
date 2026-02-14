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

def build_hole_map(h):
    N = epoch_size(h)
    amplitude = np.zeros(N, dtype=np.uint8)
    xmax = int(ceil(sqrt(250 * N / 90))) + 10
    for l, m, z in operators:
        for x in range(1, xmax + 1):
            y = 90 * x*x - l * x + m
            if y >= N: continue
            p = z + 90 * (x - 1)
            if p <= 0: continue
            idx = int(y)
            while idx < N:
                amplitude[idx] += 1
                idx += p
    return amplitude == 0  # True where holes (primes)

def to_image(arr, side):
    padded = np.pad(arr, (0, side*side - len(arr)), mode='constant')
    return padded.reshape(side, side)

def animate_epochs(h_min=10, h_max=100, h_step=5, canvas_side=512, fps=5, save=True):
    h_values = list(range(h_min, h_max + 1, h_step))
    frames = []
    for h in h_values:
        holes = build_hole_map(h)
        # Resample to fixed canvas (interpolate for smooth evolution)
        orig_side = int(ceil(sqrt(len(holes))))
        img = to_image(holes, orig_side)
        # Resize to canvas_side via nearest (preserves binary)
        x = np.linspace(0, orig_side-1, canvas_side)
        y = np.linspace(0, orig_side-1, canvas_side)
        X, Y = np.meshgrid(x, y)
        resized = np.round(np.interp(X.flatten(), range(orig_side), range(orig_side))).astype(int)
        resized_img = img[resized // orig_side, resized % orig_side].reshape(canvas_side, canvas_side)
        frames.append(resized_img)

    fig, ax = plt.subplots(figsize=(8,8))
    im = ax.imshow(frames[0], cmap='gray_r', interpolation='nearest')
    ax.axis('off')
    ax.set_title(f"Hole Evolution: h={h_values[0]} (Primes as white fringes)")

    def update(frame):
        im.set_array(frames[frame])
        ax.set_title(f"Hole Evolution: h={h_values[frame]} (Primes as white fringes)")
        return [im]

    anim = FuncAnimation(fig, update, frames=len(frames), interval=1000/fps, blit=True)
    if save:
        anim.save(f"hole_evolution_h{h_min}-{h_max}.gif", writer=PillowWriter(fps=fps))
        print(f"Saved animation as hole_evolution_h{h_min}-{h_max}.gif")
    plt.show()

# ==================== RUN IT ====================
if __name__ == "__main__":
    h_min = int(input("Enter min h (e.g. 10): ") or 10)
    h_max = int(input("Enter max h (e.g. 100): ") or 100)
    h_step = int(input("Enter step (e.g. 5): ") or 5)
    canvas_side = int(input("Canvas side (e.g. 512): ") or 512)
    fps = int(input("FPS (e.g. 5): ") or 5)
    animate_epochs(h_min, h_max, h_step, canvas_side, fps, save=True)