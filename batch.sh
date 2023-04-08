#!/bin/bash

start=$(date -d "2021-07-30 06:00:00" +%s)
end=$(date -d "2021-07-30 06:30:00" +%s)
increment=30

while [ $start -le $end ]
do
  date -d @$start +"%Y%m%d%H%M%S"
  start=$((start+increment))
  new_date=$(date -d @$start +"%Y%m%d%H%M%S")

  python pawr_to_mask_cappi.py /data/ra000007/zhaoyang/data/scale_database/scale-letkf-test-suite/obs/PAWR_Saitama/radar_${new_date}.dat

done
