#!/bin/bash

csv_file="KRIW_closest_times.csv"
destination_dir="pvol/"  # Make sure this directory exists.

# Ensure destination directory exists
mkdir -p "$destination_dir"

# Function to download a file if it exists
download_file() {
    local s3_path="$1"
    local destination="$2"  # Correctly capture the destination directory.
    # Use `aws s3 ls` to check if the file exists in S3
    if aws s3 ls "$s3_path" --no-sign-request > /dev/null 2>&1; then
        echo "Downloading $s3_path to ${destination}..."
        aws s3 cp "$s3_path" "${destination}" --no-sign-request  # Ensure ${destination} is used here.
    else
        echo "File does not exist: $s3_path"
    fi
}

export -f download_file

# Generate a list of S3 paths to download
s3_paths=()
while IFS= read -r filename; do
    # Skip empty filenames or invalid dates
    if [[ -z "$filename" || "$filename" == "s3://noaa-nexrad-level2/////" ]]; then
        continue
    fi

    radar=$(echo "$filename" | cut -c1-4)
    date=$(echo "$filename" | cut -c5-12)
    year=$(echo "$date" | cut -c1-4)
    month=$(echo "$date" | cut -c5-6)
    day=$(echo "$date" | cut -c7-8)

    # Construct the S3 path
    s3_path="s3://noaa-nexrad-level2/${year}/${month}/${day}/${radar}/${filename}"
    s3_paths+=("$s3_path")
done < <(tail -n +2 "$csv_file" | cut -d, -f8 | sed 's/_MDM$//')

# Use xargs to download files in parallel
printf "%s\n" "${s3_paths[@]}" | xargs -n 1 -P 10 -I {} bash -c 'download_file "$@"' _ {}

echo "Download process completed."
