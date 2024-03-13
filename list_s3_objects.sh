#!/bin/bash

radars=('KCYS' 'KRIW')
years=$(seq 2007 2008)
BUCKET_NAME="noaa-nexrad-level2"

generate_dates() {
    start_date=$1
    end_date=$2
    current_date=$(date -j -f "%Y-%m-%d" "$start_date" +%Y-%m-%d)

    while [[ "$current_date" < "$end_date" ]]; do
        echo "$current_date"
        current_date=$(date -j -v+1d -f "%Y-%m-%d" "$current_date" +%Y-%m-%d)
    done
}

for radar in ${radars[@]}; do
    for year in ${years[@]}; do
        
        # spring
        spring_start="${year}-03-01"
        spring_end="${year}-06-15"
        spring_dates=($(generate_dates $spring_start $spring_end))
        
        # fall period
        fall_start="${year}-08-15"
        fall_end="${year}-11-30"
        fall_dates=($(generate_dates $fall_start $fall_end))
        
        # combine
        dates=("${spring_dates[@]}" "${fall_dates[@]}")
        
        for date in "${dates[@]}"; do
            echo "Radar: $radar, Date: $date"
            s3_date=$(echo "$date" | sed 's/-/\//g')
            s3_path="s3://${BUCKET_NAME}/${s3_date}/${radar}/"            
            s3_date=$(echo "$date" | sed 's/-/\//g')
            s3_path="s3://${BUCKET_NAME}/${s3_date}/${radar}/"
            aws s3 ls --no-sign-request "$s3_path" | \
            awk '{print $4}' | grep '.gz$' | \
            parallel -j 0 -k echo | python script.py
        done
    done
done
