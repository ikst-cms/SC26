from pylab import *
import pylab as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as mticker

def plot_histogram_list_size(ax, folder_path, domain_size=100, cores=1, details='', threshold_for_labels=0):
    df = pd.read_csv(f'{folder_path}/custom_list_size', header=None)
    print(df)
    list_sizes = df[0].astype(int)

    counts, bins, patches = ax.hist(
        list_sizes,
        bins=40,
        rwidth=0.75,
        color="#03529B",
    )

    for count, patch in zip(counts, patches):
        if int(count) > threshold_for_labels:
            x_pos = patch.get_x() + patch.get_width() / 2
            y_pos = patch.get_height()
            ax.text(x_pos, y_pos + 0.01 * max(counts), f'{int(count)}', ha='center', va='bottom', fontsize=20)

    ax.set_ylim(0, max(counts) * 1.10)
    avg_list_size = int(df[0].mean())
    ax.set_xlabel("List Size", fontsize=24)
    ax.set_ylabel("Frequency", fontsize=24)
    ax.set_title(f"Domain Size: {domain_size}", fontsize=26, pad=8)
    ax.tick_params(axis='x', labelsize=24)
    ax.tick_params(axis='y', labelsize=24)
    print(f"done with domain size {domain_size} cores {cores}")
    return avg_list_size


details = 'quaternary'
cores = 64
domain_sizes = [1000, 10000, 100000, 1000000]
threshold_for_labels=[0,250,250,250]

fig = plt.figure(figsize=(20, 18))
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.28)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

for ax, ds, th in zip([ax1, ax2, ax3, ax4], domain_sizes, threshold_for_labels):
   
    plot_histogram_list_size(
        ax,
        folder_path=f'domain_size_scaling/{ds}',
        domain_size=ds,
        cores=cores,
        details=details,
        threshold_for_labels=th
    )

fig.savefig(f'fig-9-{details}-domain-size-histograms-subplots.png', dpi=300, bbox_inches='tight', pad_inches=0.03)
plt.show()
plt.close(fig)
