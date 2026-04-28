from pylab import *
import numpy as np
import pandas as pd
import matplotlib.ticker as mticker


def read_csv_and_return_avg(folder_path, domain_size=100, cores=1, details=''):
    df = pd.read_csv(f'{folder_path}/list_profile', header=None)
    avg_list_size = int(df[0].mean())
    avg_list_exec_time = df[1].mean()
    avg_list_combine_time = df[2].mean()
    return avg_list_size, avg_list_exec_time, avg_list_combine_time


if __name__=="__main__":
    
    details = 'quaternary'
    data = []
    for i in [1000,10000,100000]:
        avg_list_size,avg_list_exec_time,avg_list_combine_time= read_csv_and_plot(folder_path=f'domain_size_scaling/{i}', domain_size=i, cores=64, details=details)
        data.append([i,64,avg_list_size, avg_list_exec_time, avg_list_combine_time])


    df_domain_summary = pd.DataFrame(data, columns=["Domain Size", "Cores Count", "List Size", "Exec Time", "Combine Time"])

    data = []
    for i in [1,2,4,8,16,32,64]:
        avg_list_size, avg_list_exec_time, avg_list_combine_time = read_csv_and_plot(folder_path=f'cores_scaling/core-{i}', domain_size=10000, cores=i, details=details)
        data.append([10000,i,avg_list_size, avg_list_exec_time, avg_list_combine_time])


    df_core_summary = pd.DataFrame(data, columns=["Domain Size", "Cores Count", "List Size", "Exec Time", "Combine Time"])

    fig, ax_domain = plt.subplots(figsize=(11.2, 9.6))
    fig.subplots_adjust(top=0.84, right=0.90, left=0.11, bottom=0.13)

    d_exec, = ax_domain.plot(
        df_domain_summary["Domain Size"],
        df_domain_summary["Exec Time"],
        marker='o',
        markersize=8,
        linewidth=2.2,
        color='#4b0082',
        label='Domain Scaling - Average Parallel Execution Time per List'
    )
    d_comm, = ax_domain.plot(
        df_domain_summary["Domain Size"],
        df_domain_summary["Combine Time"],
        marker='o',
        markersize=8,
        linewidth=2.2,
        color='#8a2be2',
        label='Domain Scaling - Average Communication Time per List'
    )

    ax_domain.set_xlabel("Domain Size", fontsize=15, color='black', labelpad=8)
    ax_domain.set_ylabel("Time for Domain Scaling (sec)", fontsize=15, color='black')
    ax_domain.set_xscale('log', base=2)
    ax_domain.set_xticks(sorted(df_domain_summary["Domain Size"].unique()))
    ax_domain.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
    ax_domain.tick_params(axis='x', colors='black', labelsize=14)
    ax_domain.tick_params(axis='y', colors='black', labelsize=14)

    # Overlay second axis sharing the same drawing area (top/right only)
    ax_core = fig.add_axes(ax_domain.get_position(), frameon=False)
    ax_core.patch.set_alpha(0.0)
    ax_core.xaxis.set_label_position('top')
    ax_core.xaxis.tick_top()
    ax_core.yaxis.set_label_position('right')
    ax_core.yaxis.tick_right()
    ax_core.spines['left'].set_visible(False)
    ax_core.spines['bottom'].set_visible(False)
    ax_core.spines['top'].set_color('black')
    ax_core.spines['right'].set_color('black')

    c_exec, = ax_core.plot(
        df_core_summary["Cores Count"],
        df_core_summary["Exec Time"],
        marker='^',
        linestyle=':',
        markersize=8,
        linewidth=2.2,
        color='#0b8f78',
        label='Core Scaling - Average Parallel Execution Time per List'
    )
    c_comm, = ax_core.plot(
        df_core_summary["Cores Count"],
        df_core_summary["Combine Time"],
        marker='^',
        linestyle=':',
        markersize=8,
        linewidth=2.2,
        color='#2fbf9f',
        label='Core Scaling - Average Communication Time per List'
    )

    ax_core.set_xlabel("Number of Cores", fontsize=15, color='black', labelpad=8)
    ax_core.set_ylabel("Time for Core Scaling (sec)", fontsize=15, color='black')
    ax_core.set_xscale('log', base=2)
    ax_core.set_xticks(sorted(df_core_summary["Cores Count"].unique()))
    ax_core.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
    ax_core.tick_params(axis='x', colors='black', labelsize=14, pad=3)
    ax_core.tick_params(axis='y', colors='black', labelsize=14)
    ax_core.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    ax_domain.minorticks_on()
    ax_domain.grid(which='major', linestyle='--', color='#bfbfbf', linewidth=0.9, alpha=0.8)
    ax_domain.grid(which='minor', linestyle='--', color='#bfbfbf', linewidth=0.6, alpha=0.45)

    fig.suptitle("Parallel Execution Time and Communication Time\nSystem: IKST_64 | Model: Quaternary", fontsize=19, y=0.98)

    all_handles = [d_exec, d_comm, c_exec, c_comm]
    ax_domain.legend(
        all_handles,
        [h.get_label() for h in all_handles],
        loc='upper center',
        frameon=True,
        facecolor='white',
        edgecolor='gray',
        shadow=True,
        fontsize=12
    )

    plt.savefig(f"fig-7-{details}_combined_domain_and_core_scaling_dual_axes.png", dpi=300, bbox_inches='tight', pad_inches=0.03)
    plt.show()
    plt.close(fig)
