#!/bin/bash

# Ensure the script stops on the first error
set -e

# Check for config file argument
if [ "$#" -ne 1 ]; then
    echo "Usage: sh markerpen_workflow.sh <config_file>"
    exit 1
fi

CONFIG_FILE=$1

# Run the R script
echo "Running the R script with config file: $CONFIG_FILE"
Rscript markerpen_run.R $CONFIG_FILE

# Run the Python script
echo "Running the Python script with config file: $CONFIG_FILE"
python post_marker_pen.py --config $CONFIG_FILE

echo "Workflow completed successfully!"