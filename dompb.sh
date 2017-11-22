#---------------------------shell file for running mpb and extracting data for density-of-states calculation---------------------------

# This is a script for DOS-calculation programs. For more information, please refer to our website:
# https://github.com/boyuanliuoptics/DOS-calculation/edit/master/DOS_GGR.m

#---------------------------------------
#  set file names -- user changing part
#---------------------------------------

# file names setting
file_band_ctl="3Dexample.ctl"	# 3D example
# file_band_ctl="2Dexample.ctl"	# 2D example
data_bandmap="bandmap.data"
data_bandline="bandline.data"
file_output_line="band.txt"
file_output_map_frequency="frequency_GGR.txt"
file_output_map_velocity="velocity_GGR.txt"
file_output_map_kfrequency="frequency_Tr.txt"

#---------------------------------------
#  run mpb #1 -- frequency band in line
#---------------------------------------

# if you want only want to calculate DOS and skip #1, change 'true' to any other word in '*.ctl' file.
Theswitch=$(grep 'switch->' $file_band_ctl | cut -d '>' -f 2- | sed 's/)//g')
if [ $Theswitch = true ]; then
	# compute frequency band on lines for drawing band structure 
	mpb Zone?=false $file_band_ctl > $data_bandline 

	#  data extraction and process #1
	grep freqs $data_bandline | tail -n +2 | cut -d ',' -f 3-5,7- | cut -c 2- | sed 's/,\ /\ /g' | sed 's/$/\r/' > $file_output_line
fi

#-----------------------------------------
#  run mpb #2 -- frequency band in volumn
#-----------------------------------------

# compute frequency band in the whole Brillouin zone for DOS calculation
mpb Zone?=true $file_band_ctl > $data_bandmap

#---------------------------------
#  data extraction and process #2
#---------------------------------

Thechoice=$(grep 'define-param GGR-Tr?' $file_band_ctl | cut -d ' ' -f 3- | sed 's/)//g')
if [ $Thechoice = true ]; then
	grep freqs $data_bandmap | tail -n +2 | cut -d ',' -f 3-5,7- | cut -c 2- | sed 's/,\ /\ /g' | sed 's/$/\r/' > $file_output_map_kfrequency
else
	grep freqs $data_bandmap | tail -n +2 | cut -d ',' -f 7- | cut -c 2- | sed 's/,\ /\n/g' | sed 's/$/\r/' > $file_output_map_frequency
	grep velocity $data_bandmap | cut -d ',' -f 3- | cut -c 4- | sed 's/),\ #(/\n/g' |tr -d ')' | sed 's/$/\r/' > $file_output_map_velocity
fi

