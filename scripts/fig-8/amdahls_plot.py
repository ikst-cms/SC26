"""
Two-Parameter Amdahl's Law Speedup Plot Generator
Reads from extracted_mc_data.csv and generates fits accounting 
for serial fraction (ps) and communication overhead (pc).
Speedup is based on Time per Metropolis Step.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Headless backend to prevent display errors
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# ── 1. Configuration & Setup ──────────────────────────────────────────────────
CSV_FILE = 'mc_scaling_summary.csv'

# Read the CSV
df = pd.read_csv(CSV_FILE)

# Create an output directory for the 2-parameter plots
OUT_DIR = 'scaling_plots_amdahl_2param'
os.makedirs(OUT_DIR, exist_ok=True)

# ── 2. Two-Parameter Amdahl Fit Function ──────────────────────────────────────
def amdahl_2param(N, ps, pc):
    """
    Two-Parameter Amdahl's Law
    ps: serial fraction
    pc: communication overhead penalty
    N: number of processors/cores
    """
    # Prevent division by zero if N is an array containing 0 (safeguard)
    return 1.0 / (ps + (1.0 - ps) / N + pc * N)

# ── 3. Data Preparation ───────────────────────────────────────────────────────
# Read core count from the newly added column (handles 'cores' or 'core')
core_col = 'cores' if 'cores' in df.columns else 'core'
df['core'] = df[core_col].astype(float)

# Extract system prefix from folder name (e.g., '80-80-80' from '80-80-80-16_core')
if 'core_folder' in df.columns:
    df['system'] = df['core_folder'].str.extract(r'(.*?)-\d+_core')[0]
else:
    df['system'] = 'Default_System'

# Calculate Time per Metropolis Step (if not already present)
if 'time_per_step' not in df.columns:
    df['time_per_step'] = df['metropolis_time_hrs'] / df['metropolis_steps']

# Drop missing values and invalid times
df = df.dropna(subset=['core', 'system'])
df = df[df['time_per_step'] > 0]
df = df.sort_values('core').reset_index(drop=True)

# ── 4. Fit and Plot ───────────────────────────────────────────────────────────
systems = df['system'].unique()

for sys in systems:
    sys_df = df[df['system'] == sys].copy()
    
    # Find baseline time (1 core) for this specific system
    t1_rows = sys_df[sys_df['core'] == 1]
    if t1_rows.empty:
        print(f"Error: No 1-core baseline found for {sys}. Cannot compute speedup.")
        continue

    t1 = t1_rows['time_per_step'].mean()

    # Calculate speedup: S(N) = T(1) / T(N) using time per step
    sys_df['speedup'] = t1 / sys_df['time_per_step']

    cores = sys_df['core'].values
    speedup = sys_df['speedup'].values

    # Fit Two-Parameter Amdahl's Law
    try:
        # p0 provides initial guesses: 1% serial, tiny comm overhead
        popt, _ = curve_fit(amdahl_2param, cores, speedup, 
                            p0=[0.01, 1e-5], bounds=([0, 0], [1, 1]), maxfev=5000)
        ps, pc = popt
        
        # ── Plotting ──
        fig, ax = plt.subplots(figsize=(8, 8), dpi=120)

        # Generate a smooth array of x-values for plotting lines
        c_fit = np.linspace(1, cores.max(), 500)

        # Ideal scaling: Black Dashed Line
        ax.plot(c_fit, c_fit, 
                color='black', linestyle='--', linewidth=1.5, 
                label='Ideal scaling', zorder=1)

        # Amdahl fit: Black Solid Line
        ax.plot(c_fit, amdahl_2param(c_fit, ps, pc),
                color='black', linestyle='-', linewidth=2.0,
                label='Amdahl 2-parameter fit', zorder=2)

        # Measured data: Black Dots (linestyle='')
        ax.plot(cores, speedup,
                color='black', marker='o', markersize=7,
                linewidth=0, linestyle='',  
                label='Measured speedup', zorder=3)

        # ── Formatting & Node Boundaries ──
        # Add vertical dotted lines for every 52-core node
        max_nodes = int(cores.max() // 52) + 1
        for n_nodes in range(1, max_nodes + 1):
            nb = 52 * n_nodes
            if nb > cores.max():
                break
            ax.axvline(nb, color='0.72', linestyle=':', linewidth=1.2, zorder=0)

        # Labels
        ax.set_xlabel('Number of Cores  $N$', fontsize=12)
        ax.set_ylabel('Speedup  $S(N)$', fontsize=12)
        
        # Title with ps as percentage, and dynamic pc scaling
        ps_pct = ps * 100
        
        if pc < 1e-15:  # Handle 0 or effectively zero due to float precision
            pc_str = "$p_c$ = 0"
        else:
            exponent = int(np.floor(np.log10(pc)))
            pc_coeff = pc / (10**exponent)
            pc_str = f"$p_c$ = {pc_coeff:.4f} $\\times 10^{{{exponent}}}$"
        
        title_str = (f"Speedup Scaling: {sys} (Time per Metropolis Step)\n"
                     f"$p_s$ = {ps_pct:.4f}%   |   {pc_str}")
        ax.set_title(title_str, fontsize=14, pad=15)

        # Limits and Grid - True square diagonal mapping
        max_limit = max(cores.max(), speedup.max()) + 20
        ax.set_xlim(0, max_limit)
        ax.set_ylim(0, max_limit)

        ax.grid(True, linestyle='--', alpha=0.5)

        # Legend placement
        ax.legend(loc='upper left', frameon=True, fontsize=10, shadow=True)

        # ── Save and Close ──
        plt.tight_layout()
        filename = f"{OUT_DIR}/amdahl2_{sys}_time_per_step.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close(fig) # Close the figure to free up memory
        
        print(f"Successfully generated plot: {filename}")

    except Exception as e:
        print(f"Curve fit failed for {sys}: {e}")