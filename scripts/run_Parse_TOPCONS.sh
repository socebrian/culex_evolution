#!/bin/bash
# Run TOPCONS Python parser for all zip files
# this code is used so the outputs from the python script are sent to different folders

# Loop over each zip file
for zip in *.zip; do
    # Skip if no .zip files found
    [ -e "$zip" ] || { echo "No zip files found."; exit 1; }

    base=$(basename "$zip" .zip)
    echo "=== Running $zip ==="

    # Make a separate folder for each job
    mkdir -p "$base"

    (
        cd "$base"
        # Call your Python script, passing the zip as argument
        python ../Parse_TOPCONS.py "../$zip" 2>&1 | tee "${base}_$(date +%Y%m%d_%H%M).log"
    )

    echo "=== Done $zip ==="
    echo
done