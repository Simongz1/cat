from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation
import numpy as np

# Read the output relative to this script, regardless of the current directory.
csv_path = Path(__file__).with_name("positions.csv")

# Animation settings.
output_path = Path(__file__).with_name("particles.mp4")
fps = 10
frame_stride = 25  # Set to 1 to include every recorded timestep.

def particles_per_timestep(path):
    """Count the first timestep without loading the complete CSV."""
    with path.open() as csv_file:
        next(csv_file)
        first_time = next(csv_file).split(",", 1)[0]
        count = 1
        for line in csv_file:
            if line.split(",", 1)[0] != first_time:
                return count
            count += 1
    return count


def load_plot_data(path, rows_per_timestep, stride):
    """Load animation frames only; retain one KE value per timestep."""
    frames = []
    energy_times = []
    energies = []
    selected_rows = []
    selected_time = None

    with path.open() as csv_file:
        next(csv_file)
        for row_index, line in enumerate(csv_file):
            particle_index = row_index % rows_per_timestep
            timestep_index = row_index // rows_per_timestep

            # KE is repeated for every particle, so parse it only once.
            if particle_index == 0:
                columns = line.rstrip().split(",")
                energy_times.append(float(columns[0]))
                energies.append(float(columns[4]))

            if timestep_index % stride != 0:
                continue

            columns = line.rstrip().split(",")
            if particle_index == 0:
                selected_time = float(columns[0])
                selected_rows = []
            selected_rows.append(
                (float(columns[1]), float(columns[2]), float(columns[3]))
            )

            if particle_index == rows_per_timestep - 1:
                frame = np.asarray(selected_rows, dtype=np.float32)
                frames.append((selected_time, frame))

    return frames, np.asarray(energy_times), np.asarray(energies)


rows_per_timestep = particles_per_timestep(csv_path)
frames, energy_times, energies = load_plot_data(
    csv_path, rows_per_timestep, frame_stride
)
if not frames:
    raise ValueError(f"No particle data found in {csv_path}")

frame_times = np.asarray([time for time, _ in frames])

particle_mass = 1000e-6
boltzmann_constant = 1.380649e-23


def temperature_from_speed(speed):
    kinetic_energy = 0.5 * particle_mass * np.square(speed)
    return 2.0 * kinetic_energy / (3.0 * boltzmann_constant)


temperature_min = min(
    np.min(temperature_from_speed(frame[:, 2])) for _, frame in frames
)
temperature_max = max(
    np.max(temperature_from_speed(frame[:, 2])) for _, frame in frames
)
if np.isclose(temperature_min, temperature_max):
    temperature_max = temperature_min + 1.0

initial_frame = frames[0][1]

# Keep the camera locked to the initial configuration so escaped particles
# do not force the movie to zoom out.
x_min, x_max = np.min(initial_frame[:, 0]), np.max(initial_frame[:, 0])
y_min, y_max = np.min(initial_frame[:, 1]), np.max(initial_frame[:, 1])
x_padding = 0.1 * max(x_max - x_min, 1.0)
y_padding = 0.1 * max(y_max - y_min, 1.0)

initial_temperature = temperature_from_speed(initial_frame[:, 2])

fig, ax = plt.subplots()
particles = ax.scatter(
    initial_frame[:, 0],
    initial_frame[:, 1],
    c=initial_temperature,
    s=10,
    cmap="viridis",
    vmin=temperature_min,
    vmax=temperature_max,
    edgecolors="none",
)

colorbar = fig.colorbar(particles, ax=ax)
colorbar.set_label("Temperature (K)")

#ax.set_xlim(x_min - x_padding, x_max + x_padding)
#ax.set_ylim(y_min - y_padding, y_max + y_padding)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal", adjustable="box")
title = ax.set_title(f"Particle positions at t = {frame_times[0]:g}")
fig.tight_layout()


def update(frame_index):
    time, frame = frames[frame_index]
    frame_temperature = temperature_from_speed(frame[:, 2])

    particles.set_offsets(frame[:, :2])
    particles.set_array(frame_temperature)
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
    bitrate=500,
)
animation.save(output_path, writer=writer, dpi=750)
plt.close(fig)

print(f"Saved {len(frame_times)} frames to {output_path}")

fig, ax = plt.subplots()
ax.plot(energy_times, energies, "k")
ax.set_xlabel("time")
ax.set_ylabel("kinetic energy")
plt.show()
