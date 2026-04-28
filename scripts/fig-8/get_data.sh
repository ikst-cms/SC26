#!/bin/bash

output_csv="mc_scaling_summary.csv"

# CSV Header updated for the new requirements
echo "folder_name,atoms,cores,metropolis_steps,metropolis_time_hrs" > "$output_csv"

# Loop through all directories starting with 60-60-60-
for dir in 60-60-60-*/; do
  # Get the clean folder name without the trailing slash
  dir_name=$(basename "$dir")
  
  # Assuming MC.log is directly inside this folder. 
  # If it is inside a 100000 subfolder like before, change this to: "$dir/100000/MC.log"
  mc_file="$dir/MC.log"
  
  if [ -f "$mc_file" ]; then
    # Extract number of atoms (active sites)
    # Looks for "Found X active sites out of" and grabs the X
    atoms=$(grep "active sites out of" "$mc_file" | head -n 1 | awk '{print $2}')
    
    # Extract number of cores (processes)
    # Looks for "Loop Parallelization on X processes" and grabs the X
    cores=$(grep "Loop Parallelization on" "$mc_file" | head -n 1 | grep -oE '[0-9]+ processes' | grep -oE '[0-9]+')
    
    # Extract and SUM all Metropolis steps
    metro_steps=$(grep "Metropolis equilibriation steps" "$mc_file" \
                  | grep -oE '[0-9]+ Metropolis' \
                  | grep -oE '[0-9]+' \
                  | awk '{sum += $1} END {if (NR > 0) print sum}')
    
    # Extract and SUM all Metropolis times
    metro_time=$(grep "Metropolis equilibriation steps" "$mc_file" \
                 | grep -oE '[0-9]+\.[0-9]+ hrs' \
                 | grep -oE '[0-9]+\.[0-9]+' \
                 | awk '{sum += $1} END {if (NR > 0) printf "%.6f\n", sum}')

    # Append all extracted values to the CSV
    echo "$dir_name,${atoms:-N/A},${cores:-N/A},${metro_steps:-N/A},${metro_time:-N/A}" >> "$output_csv"
  else
    echo "Warning: MC.log not found in $dir_name"
  fi
done

echo "Extraction complete. Results saved to $output_csv"