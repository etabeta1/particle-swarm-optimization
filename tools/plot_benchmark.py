#!/usr/bin/env python3
import json
import matplotlib.pyplot as plt
import sys

with open(sys.argv[1], 'r') as f:
    data = json.load(f).values()

max_cpus = max([b['threads'] for b in data])
particle_sizes = set(b['particles_per_type'] for b in data)

xs = range(1, max_cpus + 1)

# Draw linear plot
for ps in particle_sizes:
    ys = [s['time_seconds'] for s in sorted([b for b in data if b['particles_per_type'] == ps], key=lambda b: b['threads'])]
    plt.plot(xs, ys, label=f'Particles per type: {ps}')

plt.xlabel('Number of CPUs')
plt.ylabel('Time (seconds)')
plt.title(sys.argv[3] + " (linear scale)")
plt.legend()
plt.grid(True)
plt.savefig(sys.argv[2].replace(".png", "_linear.png"))

# Draw log-log plot
plt.clf()
for ps in particle_sizes:
    ys = [s['time_seconds'] for s in sorted([b for b in data if b['particles_per_type'] == ps], key=lambda b: b['threads'])]
    plt.loglog(xs, ys, label=f'Particles per type: {ps}')

# Draw reference log-log line
ref_xs = [1, max_cpus]
ref_ys = [ys[0], ys[0] / max_cpus]

plt.loglog(ref_xs, ref_ys, 'k--', label='Ideal Scaling')

# Draw log2-log2 line
ref_ys_log2 = [ys[0], ys[0] / (max_cpus ** 0.5)]
plt.loglog(ref_xs, ref_ys_log2, 'r--', label='Log2 Scaling')

plt.xlabel('Number of CPUs (log scale)')
plt.ylabel('Time (seconds, log scale)')
plt.title(sys.argv[3] + " (log-log scale)")
plt.legend()
plt.grid(True, which="both", ls="--")
plt.savefig(sys.argv[2].replace(".png", "_loglog.png"))
