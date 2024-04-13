#!/bin/bash

# radars=(
# KABR KABX KAKQ KAMA KAMX KAPX KARX KATX KBBX KBGM
# KBHX KBIS KBLX KBMX KBOX KBRO KBUF KBYX KCAE KCBW
# KCBX KCCX KCLE KCLX KCRP KCYS KDAX KDDC KDFX KDGX
# KDIX KDLH KDMX KDOX KDTX KDVN KDYX KEAX KEMX KENX
# KEOX KEPZ KESX KEVX KEWX KEYX KFCX KFDR KFDX KFFC
# KFSD KFSX KFTG KFWS KGGW KGJX KGLD KGRB KGRK KGRR
# KGSP KGWX KGYX KHDX KHGX KHNX KHPX KHTX KICT KICX
# KILN KIND KINX KIWA KIWX KJAX KJGX KJKL KLBB KLCH
# KLGX KLIX KLNX KLOT KLRX KLSX KLTX KLVX KLWX KLZK
# KMAF KMAX KMBX KMHX KMKX KMLB KMOB KMPX KMQT KMRX
# KMSX KMTX KMUX KMVX KNKX KNQA KOAX KOHX KOKX KOTX
# KPAH KPBZ KPDT KPOE KPUX KRAX KRGX KRIW KRLX KRTX
# KSFX KSGF KSHV KSJT KSOX KSRX KTBW KTFX KTLH KTLX
# KTWX KUDX KUEX KVAX KVBX KVNX KVTX KVWX KYUX PABC
# PACG PAEC PAHG PAIH PAKC PAPD PGUA PHKI PHKM PHMO
# PHWA RKJK RKSG RODN TJUA KCXX KILX KTYX KMXX KJAN
# KOUN KRMX
# )

#radars=($1)
#start_year=$2
#end_year=$3

radars=(${RADAR})
start_year=${START_YEAR}
end_year=${END_YEAR}

years=$(seq $start_year $end_year)

BUCKET_NAME="noaa-nexrad-level2"

generate_dates() {
    start_date=$1
    end_date=$2
    current_date=$(date -d "$start_date" +%Y-%m-%d)
    end_date_timestamp=$(date -d "$end_date" +%s)

    while : ; do
        current_date_timestamp=$(date -d "$current_date" +%s)

        # Break the loop if current_date is greater than end_date
        if [[ "$current_date_timestamp" -gt "$end_date_timestamp" ]]; then
            break
        fi

        echo "$current_date"
        # Increment the date by one day and format it.
        current_date=$(date -d "$current_date + 1 day" +%Y-%m-%d)
    done
}


for radar in ${radars[@]}; do
    for year in ${years[@]}; do
        
        # spring
        spring_start="${year}-01-01"
        spring_end="${year}-06-01"
        spring_dates=($(generate_dates $spring_start $spring_end))
        
        # fall period
        fall_start="${year}-06-02"
        fall_end="${year}-12-31"
        fall_dates=($(generate_dates $fall_start $fall_end))
        
        # combine
        #dates=("${fall_dates[@]}")
        dates=("${spring_dates[@]}" "${fall_dates[@]}")
        
        for date in "${dates[@]}"; do
            echo "Radar: $radar, Date: $date"
            s3_date=$(echo "$date" | sed 's/-/\//g')
            s3_path="s3://${BUCKET_NAME}/${s3_date}/${radar}/"

            # Use aws s3 ls to list files and capture the output
            s3_files=$(aws s3 ls --no-sign-request "$s3_path" | awk '{print $4}')

            # Check if s3_files is not empty
            if [[ -n "$s3_files" ]]; then
                # If not empty, there are files, so proceed with the script
                echo "$s3_files" | parallel -j 0 -k echo | python script.py
            else
                # If empty, skip this date for this radar
                echo "No data for Radar: $radar, Date: $date"
            fi
        done
    done
done
