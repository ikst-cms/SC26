from pylab import *
import numpy as np
import pandas as pd
import matplotlib.ticker as mticker

def plot_histogram_list_size(ax, folder_path, swap_dis=0.9, cores=1, details='', threshold_for_labels = 0):
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
    ax.set_xlabel("List Size", fontsize=22)
    ax.set_ylabel("Frequency", fontsize=22)
    ax.set_title(f"Threshold Swap Distance: {swap_dis}", fontsize=26, pad=8)
    ax.tick_params(axis='x', labelsize=22)
    ax.tick_params(axis='y', labelsize=22)
    return avg_list_size


details = 'quaternary_10000'
cores = 64
swap_dis = [0.9, 3]
threshold_for_labels=[250,500]

fig, axes = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)

for ax, ds, th in zip(axes, swap_dis, threshold_for_labels):
    plot_histogram_list_size(
        ax,
        folder_path=f'swap-dist-{ds}',
        swap_dis=ds,
        cores=cores,
        details=details,
        threshold_for_labels=th
    )

fig.savefig(
    f'fig-10-{details}-swap-dist-histograms.png',
    dpi=300,
    bbox_inches='tight',
    pad_inches=0.03
)
plt.show()
plt.close(fig)
