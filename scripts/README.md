


## Purpose
This folder contains driver scripts and plotting utilities used to reproduce the figures in the performance study. Below are the folder layout and per-figure run instructions.

## Layout
- `fig-5,6/` — domain-time scaling driver and plotting script (`domain_time_scaling.sh`, `domain_time_plots.py`)
- `fig-7,9/` — domain scaling and core-scaling drivers and plotting scripts (`sample_domain_scaling.sh`, `sample_cores_scaling.sh`, `domain_scaling_job.sh`, `core_scaling_job.sh`, `plot-fig-7.py`, `plot-fig-9.py`)
- `fig-8/` — speedup scaling driver(s) and plotting helpers
- `fig-10/` — swap-distance driver and histogram plotting (`sample_swap_dist_scaling.sh`, `plot-fig-10.py`, `swap-dist-*/` example inputs)

## How to run 

### Fig. 5 & 6 (MC Time/Step and MC Time/Flip)

- Folder: `scripts/fig-5,6`
- Purpose: Measure MC time per step and per flip across domain sizes and node configurations.
- Domain sizes: 1e3, 1e4, 1e5.
- Core configurations: single-node (1 and 8 cores), two-node (104 cores), six-node (312 cores) — sample PBS scripts encode these allocations.
- Script to run: `domain_time_scaling.sh` (driver that launches the three PBS/job variants).
- Per-run outputs: `MC.log` (contains timing), which `domain_time_scaling.sh` aggregates into `mc_summary.csv`.
- Post-process / plot: `domain_time_plots.py mc_summary.csv` produces the MC Time/Step and MC Time/Flip figures.

Commands

```bash
cd scripts/fig-5,6
./domain_time_scaling.sh
python3 domain_time_plots.py mc_summary.csv
```

### Fig. 7 & 9 (List profiling and list-size histograms)

- Folder: `scripts/fig-7,9`
- Purpose: (a) measure per-list execution and MPI communication times; (b) collect list-size distributions.
- Domain scaling: domain sizes 1e3, 1e4, 1e5, 1e6 on a 64-core node. Use `sample_domain_scaling.sh` to run these. `domain_scaling_job.sh` is a sample PBS job script file.
- Core scaling: for domain = 1e4, scale cores across 1,2,4,8,16,32,64 using `sample_cores_scaling.sh`. `core_scaling_job.sh` is a sample PBS job script file.
- It  produces `list_profile` (per-list times + comm time) and `custom_list_size` (list sizes) CSVs.
- Post-process / plot: `plot-fig-7.py` consumes `list_profile` files to compute average communication (MPI Allreduce) vs compute times; `plot-fig-9.py` consumes `custom_list_size` to build histograms.

Commands

```bash
cd scripts/fig-7,9
./sample_domain_scaling.sh    # domain scaling at 64 cores
./sample_cores_scaling.sh     # core scaling for domain 1e4
python3 plot-fig-7.py
python3 plot-fig-9.py
```

### Fig. 8 (Speedup)

- Folder: `scripts/fig-8`
- Purpose: Strong-scaling (speedup) experiment for a HEA 6-component dataset (example: 60×60×60).
- Core range: 1 up to 512 cores (scripts submit jobs across increasing core counts; use cluster PBS variants as required).
- Outputs: timing logs per run; post-process to compute speedup relative to single-core baseline.
- Plot: `amdahls.py` or other helper after collecting timings.

Commands

```bash
cd scripts/fig-8
./job.sh
python3 amdahls.py
```

### Fig. 10 (Swap distance sensitivity)


- Folder: `scripts/fig-10`
- Purpose: Compare list-size distributions for different swap-distance thresholds.
- Dataset: quaternary dataset at domain size 1e4 (example inputs in `swap-dist-0.9/` and `swap-dist-3/`).
- Required cores: at least 32 cores (script assumes mpiexec -np 32 in its driver).
- Script: `sample_swap_dist_scaling.sh` runs CPyMonC in each `swap-dist-*` folder. The run logs list sizes in `custom_list_size` files.
- Plot: `plot-fig-10.py` reads each `swap-dist-*/custom_list_size` and generates comparative histograms.

Commands
```bash
cd scripts/fig-10
./sample_swap_dist_scaling.sh
python3 plot-fig-10.py
```




