#!/bin/bash

state="${1:-0}"

if [[ "$state" != "0" && "$state" != "1" ]]; then
    echo "Usage: $0 [0|1]"
    echo "  0 -> build GIF from state0 snapshots"
    echo "  1 -> build GIF from state1 snapshots"
    exit 1
fi

for file in output/wf_*.dat
do
    echo "Plotting $file"
    gnuplot -c plot_density.gscript "$file"
done

mapfile -t frames < <(printf '%s\n' output/wf_*_state"${state}".png | sort -V)

if [[ ${#frames[@]} -eq 0 || ! -e "${frames[0]}" ]]; then
    echo "No PNG frames found for state${state}"
    exit 1
fi

echo "Creating GIF for state${state}"
magick -delay 30 -loop 0 "${frames[@]}" "output/wf_state${state}.gif"
