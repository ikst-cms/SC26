"""
Consolidated Scaling Plots Generator
Reads from mc_summary.csv in the current directory.
Generates both Time per Flip and Domain Size scaling plots as log-log plots.

Usage: python generate_scaling_plots.py
"""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # Headless backend
import matplotlib.pyplot as plt
import os

# ── 1. Configuration & Setup ──────────────────────────────────────────────────
CSV_FILE = "mc_summary.csv"

# Check if mc_summary.csv exists in current directory
if not os.path.exists(CSV_FILE):
    print(f"Error: {CSV_FILE} not found in current directory")
    exit(1)

print(f"Processing {CSV_FILE}...")

# Read the CSV
df = pd.read_csv(CSV_FILE)

# Ignore supercell size of 100
df = df[df["supercell_size"] != 100]

# Calculate Number of Atoms (BCC Lattice = 2 * supercell_size)
df["atoms"] = df["supercell_size"] * 2

# Calculate MC Time/Step (sec)
df["mc_time_per_step_sec"] = (df["total_metro_time_hrs"] * 3600) / df[
    "total_metro_steps"
]

# time_per_flip_secs is already in the CSV, but ensure it's numeric
df["time_per_flip_sec"] = pd.to_numeric(df["time_per_flip_secs"], errors="coerce")

# Create output directories for plots
out_dir = "scaling_plots"
os.makedirs(out_dir, exist_ok=True)

# Get unique cores (representing different MPI task counts)
cores = sorted(df["model"].dropna().unique())

# Extract numeric part from core names (e.g., 'quaternary-core-1' -> 1)
core_numbers = {}
for core in cores:
    try:
        core_num = int(core.split('-')[-1])
        core_numbers[core] = core_num
    except:
        core_numbers[core] = core

# Sort cores by their numeric values
cores = sorted(cores, key=lambda x: core_numbers[x])

# ── 2. Generate Plots ─────────────────────────────────────────────────────

# Select representative cores (max 4)
if len(cores) > 4:
    indices = np.linspace(0, len(cores) - 1, 4).astype(int)
    plot_cores = [cores[i] for i in indices]
else:
    plot_cores = cores

# Define aesthetics for the different nodes using a colormap
cmap = plt.cm.viridis
colors = cmap(np.linspace(0, 0.9, len(plot_cores)))
markers = ["o", "s", "^", "D", "v", "<", ">", "p", "*", "h", "H", "+", "x"] * 5

# ── 2.1 Domain Size Scaling Plot ──────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 8), dpi=120)

for idx, core in enumerate(plot_cores):
    core_df = df[df["model"] == core].sort_values("atoms")
    if core_df.empty:
        continue

    # Filter out any N/A or invalid values
    core_df = core_df[core_df["mc_time_per_step_sec"].notna()]
    core_df = core_df[core_df["mc_time_per_step_sec"] > 0]

    if core_df.empty:
        continue

    ax.plot(
        core_df["atoms"],
        core_df["mc_time_per_step_sec"],
        marker=markers[idx],
        color=colors[idx],
        linestyle="-",
        linewidth=2,
        markersize=8,
        label=f"{core_numbers[core]} MPI Task(s)",
    )

# Formatting
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Number of Atoms", fontsize=12)
ax.set_ylabel("MC Time/Step (sec)", fontsize=12)
ax.set_title("Domain Size Scaling", fontsize=14, pad=15)
ax.grid(True, which="both", linestyle="--", alpha=0.5)
ax.legend(loc="upper left", frameon=True, fontsize=10, shadow=True)

plt.tight_layout()
domain_filename = os.path.join(out_dir, "domain_scaling.png")
plt.savefig(domain_filename, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {domain_filename}")

# ── 2.2 Time per Flip Scaling Plot ────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 8), dpi=120)

for idx, core in enumerate(plot_cores):
    core_df = df[df["model"] == core].sort_values("atoms")
    if core_df.empty:
        continue

    # Filter out any N/A or invalid values
    core_df = core_df[core_df["time_per_flip_sec"].notna()]
    core_df = core_df[core_df["time_per_flip_sec"] > 0]

    if core_df.empty:
        continue

    ax.plot(
        core_df["atoms"],
        core_df["time_per_flip_sec"],
        marker=markers[idx],
        color=colors[idx],
        linestyle="-",
        linewidth=2,
        markersize=8,
        label=f"{core_numbers[core]} MPI Task(s)",
    )

# Formatting
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Number of Atoms", fontsize=12)
ax.set_ylabel("MC Time / Flip (seconds)", fontsize=12)
ax.set_title("Time per Flip Scaling", fontsize=14, pad=15)
ax.grid(True, which="both", linestyle="--", alpha=0.5)
ax.legend(loc="best", frameon=True, fontsize=10, shadow=True)

plt.tight_layout()
flip_filename = os.path.join(out_dir, "flip_scaling.png")
plt.savefig(flip_filename, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {flip_filename}")

print("\nSuccessfully generated the plots.")
