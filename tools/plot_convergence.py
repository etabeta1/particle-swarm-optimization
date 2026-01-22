#/usr/bin/env python3
from collections import defaultdict
import matplotlib.pyplot as plt
import sys

# Read the file
x_values = defaultdict(list)

with open(sys.argv[1], 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 5:
            y = float(parts[2])
            x = float(parts[3])
            type_val = int(parts[4])
            
            if type_val == 0 and x < 256:
                x_values[x].append(y)

x_sorted = sorted(x_values.keys())

minimums = [min(x_values[x]) for x in x_sorted]
cumulative_minumums = []
current_min = float('inf')
for val in minimums:
    if val < current_min:
        current_min = val
    cumulative_minumums.append(current_min)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(x_sorted, cumulative_minumums, label='Global best')
plt.plot(x_sorted, minimums, label='Minimum of local best')
plt.title('Convergence Plot of Objective Function')
plt.xlabel('Iteration')
plt.ylabel('Objective function value')
plt.legend()
plt.grid(True)
plt.savefig(sys.argv[2])