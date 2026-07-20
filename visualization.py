from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation
import numpy as np

# Read the output relative to this script, regardless of the current directory.
csv_path = Path(__file__).with_name("positions.csv")
data = np.genfromtxt(csv_path, delimiter=",", names=True)

# Animation settings.
output_path = Path(__file__).with_name("particles.mp4")
fps = 30
frame_stride = 5  # Set to 1 to include every recorded timestep.

all_times = data["time"]
frame_times = np.unique(all_times)[::frame_stride]

speed_min = np.min(data["vmag"])
speed_max = np.max(data["vmag"])
if np.isclose(speed_min, speed_max):
    speed_max = speed_min + 1.0


def rows_at_time(time):
    start = np.searchsorted(all_times, time, side="left")
    stop = np.searchsorted(all_times, time, side="right")
    return data[start:stop]


initial_frame = rows_at_time(frame_times[0])

# Keep the camera locked to the initial configuration so escaped particles
# do not force the movie to zoom out.
x_min, x_max = np.min(initial_frame["x"]), np.max(initial_frame["x"])
y_min, y_max = np.min(initial_frame["y"]), np.max(initial_frame["y"])
x_padding = 0.1 * max(x_max - x_min, 1.0)
y_padding = 0.1 * max(y_max - y_min, 1.0)

fig, ax = plt.subplots()
particles = ax.scatter(
    initial_frame["x"],
    initial_frame["y"],
    c=initial_frame["vmag"],
    s=40,
    cmap="viridis",
    vmin=speed_min,
    vmax=speed_max,
    edgecolors="none",
)

colorbar = fig.colorbar(particles, ax=ax)
colorbar.set_label("Velocity magnitude")

#ax.set_xlim(x_min - x_padding, x_max + x_padding)
#ax.set_ylim(y_min - y_padding, y_max + y_padding)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal", adjustable="box")
title = ax.set_title(f"Particle positions at t = {frame_times[0]:g}")
fig.tight_layout()


def update(frame_index):
    time = frame_times[frame_index]
    frame = rows_at_time(time)

    particles.set_offsets(np.column_stack((frame["x"], frame["y"])))
    particles.set_array(frame["vmag"])
    title.set_text(f"Particle positions at t = {time:g}")

    return particles, title


animation = FuncAnimation(
    fig,
    update,
    frames=len(frame_times),
    interval = 1000 / fps,
    blit=True,
)

writer = FFMpegWriter(
    fps=fps,
    metadata={"title": "Particle hydrodynamics"},
    bitrate=2400,
)
animation.save(output_path, writer=writer, dpi=150)
plt.close(fig)

print(f"Saved {len(frame_times)} frames to {output_path}")

fig, ax = plt.subplots()
ax.plot(data['time'], data['KE'], 'k')
plt.show()
