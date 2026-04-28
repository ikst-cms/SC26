#!/bin/bash

output_csv="extracted_mc_data.csv"

# CSV Header simplified since there are no model or supercell folders here
echo "core_folder,nn_list_time_secs,energy_diff_time_secs,metropolis_steps,metropolis_time_hrs,flips" > "$output_csv"

# Loop through all directories in the current folder (e.g., 1, 2, 4, 8, 16...)
for core_dir in */; do
  # Remove the trailing slash to get just the folder name
  core_num=$(basename "$core_dir")
  
  # Paths to the log and list size files
  mc_file="${core_dir}MC.log"
  flips_file="${core_dir}custom_list_size"
  
  if [ -f "$mc_file" ]; then
    # Extract Nearest Neighbor list time
    nn_time=$(grep "Successfully created nearest-neighbor list" "$mc_file" | grep -oE '[0-9]+\.[0-9]+')
    
    # Extract Energy difference time
    energy_time=$(grep "Using Local calculation to find energy difference" "$mc_file" | grep -oE '[0-9]+\.[0-9]+')

    # Extract and SUM all Metropolis steps using awk
    metro_steps=$(grep "Metropolis equilibriation steps" "$mc_file" \
                  | grep -oE '[0-9]+ Metropolis' \
                  | grep -oE '[0-9]+' \
                  | awk '{sum += $1} END {if (NR > 0) print sum}')
    
    # Extract and SUM all Metropolis times (floating point addition) using awk
    metro_time=$(grep "Metropolis equilibriation steps" "$mc_file" \
                 | grep -oE '[0-9]+\.[0-9]+ hrs' \
                 | grep -oE '[0-9]+\.[0-9]+' \
                 | awk '{sum += $1} END {if (NR > 0) printf "%.6f\n", sum}')

    # Extract and SUM flips from custom_list_size
    if [ -f "$flips_file" ]; then
      total_flips=$(awk 'NF {for(i=1;i<=NF;i++) sum+=$i; has_data=1} END {if(has_data) print sum; else print "N/A"}' "$flips_file")
    else
      total_flips="N/A"
    fi

    # Append all extracted values to the CSV
    echo "$core_num,${nn_time:-N/A},${energy_time:-N/A},${metro_steps:-N/A},${metro_time:-N/A},${total_flips}" >> "$output_csv"
  fi
done

echo "Extraction complete! Results saved to $output_csv"