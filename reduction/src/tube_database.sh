#/bin/bash

#

# On current directory there must be two filterbank files, one aiming to Hydra (on.fil) an the other offset (off.fil)

#

numfilfiles=(*.fil)
numfilfiles=${#numfilfiles[@]}
if [ $numfilfiles -ne 2 ]; then
	echo 'ERROR: directory must contain two observations, one aiming to Hydra (on.fil) an the other offset (off.fil).'
	echo 'E.g.: ds256_pulcal_020171023_104155_on.fil, ds256_pulcal_020171023_104155_off.fil'
	return 1
fi

if [ -d "reduc*" ]; then
	read -p 'Existing reduction in directory, delete and make a new one? [y/n] ' -e -i 'y' delete_ans
	dmsearch="$(echo ${dmsearch} | tr 'A-Z' 'a-z')"
	if [ $delete_ans == "y" ] || [ $dmsearch == "yes" ]; then
		rm -r reduc*
		else if [ $delete_ans == "n" ] || [ $delete_ans == "no" ]; then
			return 1
		else
			echo "ERROR: existing reductions?"
			return 1
	fi
fi

# Reading observation files:
on_fil=$(ls *on.fil)
off_fil=$(ls *off.fil)

# Timestamp:
TIMESTAMP=$(readfile $onfil | grep 'Obs' | awk '{print $5}')

# Tube parameters:
period=3.0
dm=0.0

# Mask parameters:
timemask=1.0

# Folding parameters:
nsub=32
nbins=512
phase=0.0

# Creating reduction directories:
mkdir reduc_on
mkdir reduc_off

####################################
# REDUCTION STARTS, ON HYDRA:
####################################
cd reduc_on
# Making mask:
rfifind -time $timemask -zerodm -o mask_pul_cal_on ../$pulcalfil_on
# Folding:
prepfold -nsub $nsub -n $nbins -p $period -nopsearch -dm $dm -nodmsearch -phs $phase -o prepfold_$TIMESTAMP_tube_on -noxwin ./$on_fil $mask_cal_on
# Editing tube_on.pfd:
psredit -c rcvr:name=AI -c be:name=EB120 -c coord=09:18:05.651-12:05:43.99 -c type=FluxCal-On -m prepfold_$TIMESTAMP_tube_on.pfd
cp prepfold_$TIMESTAMP_tube_on.pfd ../../fcal_files/$TIMESTAMP_on.fcal


####################################
# REDUCTION STARTS, OFF HYDRA:
####################################
cd ..
cd reduc_off
# Making mask:
rfifind -time $timemask -zerodm -o mask_pul_cal_of ../$pulcalfil_of
# Folding:
prepfold -nsub $nsub -n $nbins -p $period -nopsearch -dm $dm -nodmsearch -phs $phase -o prepfold_$TIMESTAMP_tube_off -noxwin ./$off_fil $mask_cal_off
# Editing tube_off.pfd:
psredit -c rcvr:name=AI -c be:name=EB120 -c coord=09:18:05.651-14:05:43.99 -c type=FluxCal-Off -m prepfold_$TIMESTAMP_tube_on.pfd
cp prepfold_$TIMESTAMP_tube_off.pfd ../../fcal_files/$TIMESTAMP_off.fcal

