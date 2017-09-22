file_band_ctl="DG.ctl"
data_bandmap="bandmap.data"
data_bandline="bandline.data"

file_output_line="band.dat"
file_output_map_frequency="frequency_GR.txt"
file_output_map_velocity="velocity_GR.txt"
# compute frequency band on high symmetry lines
mpb Zone?=false $file_band_ctl > $data_bandline 
#mpb-split 8 Zone?=false $file_band_ctl > $data_bandline  # for multithreading, '8' is the number of threads
grep freqs $data_bandline > $file_output_line

# compute frequency band in the whole Brillouin zone
mpb Zone?=true $file_band_ctl > $data_bandmap
#mpb-split 8 Zone?=true $file_band_ctl > $data_bandmap    # for multithreading, '8' is the number of threads
grep freqs $data_bandmap | tail -n +2 | cut -d ',' -f 7- | cut -c 2- | sed 's/,\ /\n/g' | sed -e 's/$/\r/' > $file_output_map_frequency
grep velocity $data_bandmap | cut -d ',' -f 3- | cut -c 4-  | sed 's/),\ #(/\n/g' |tr -d ')' | sed -e 's/$/\r/' > $file_output_map_velocity


