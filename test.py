import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

# Create sample data
np.random.seed(42)  # For reproducibility
x = np.random.rand(20)
y = np.random.rand(20)

# Create a figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# First plot: 5x5 cm with "1cm" markers (relative to the plot)
# In a 5x5 cm plot with data range 0-1, 1cm = 0.2 data units
marker_size_in_data_units = 0.2

# For the first plot
s1 = 400  # This gives approximately the desired size in the first plot

# Create the first scatter plot
ax1.scatter(x, y, s=s1, alpha=0.7, color='blue')
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)
ax1.set_title("5x5 cm plot\nMarkers ~1cm wide")
ax1.grid(True, linestyle='--', alpha=0.7)

# Add a reference circle to visualize the "1cm" size
circle = patches.Circle((0.15, 0.85), 0.1, fill=False, edgecolor='red', linestyle='--')
ax1.add_patch(circle)
ax1.text(0.15, 0.75, '1cm equivalent\n(0.2 data units)', ha='center', fontsize=8)

# Second plot: 10x10 cm with same relative marker size
# The same 0.2 data units should still appear as "1cm" relative to plot
s2 = s1  # Same marker size parameter maintains the relative size

# Create the second scatter plot
ax2.scatter(x, y, s=s2, alpha=0.7, color='blue')
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
ax2.set_title("10x10 cm plot\nSame relative marker size")
ax2.grid(True, linestyle='--', alpha=0.7)

# Add a reference circle with the same data unit size
circle = patches.Circle((0.15, 0.85), 0.1, fill=False, edgecolor='red', linestyle='--')
ax2.add_patch(circle)
ax2.text(0.15, 0.75, '1cm equivalent\n(0.2 data units)', ha='center', fontsize=8)

# Add an explanatory note
plt.figtext(0.5, 0.01, 
            "Note: The red circles represent the same data units (0.2) in both plots.\n"
            "When the plot area expands from 5x5cm to 10x10cm, the markers maintain\n"
            "the same visual size relative to the data range (not physical size).",
            ha='center', fontsize=9)

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
plt.savefig('scatter_plot_comparison.png', dpi=100)
plt.show()