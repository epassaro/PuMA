#/bin/bash

#

# On current directory there must be two filterbank files, one aiming to Hydra (on.fil) an the second to the Celestial South Pole (off.fil)

#

numfilfiles=(*.fil)
numfilfiles=${#numfilfiles[@]}
if [ $numfilfiles -gt 2 ]; then
	echo 'ERROR: too many .fil files in this directory, must contain only pulsar and calibration observations.'
	echo 'E.g.: ds256_B1240-64_020171023_104155.fil, ds256_pul_cal_020171023_104155.fil'
	exit
fi

# Timestamp:
TIMESTAMP=$(date +"%y%m%d_%H%M%S")

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
mkdir reduc_on_$TIMESTAMP
mkdir reduc_off_$TIMESTAMP


####################################
# REDUCTION STARTS, ON HYDRA:
####################################

cd reduc_on_$TIMESTAMP

# Making mask:
rfifind -time $timemask -zerodm -o mask_pul_cal_on ../$pulcalfil_on
# Folding:
prepfold -nsub $nsub -n $nbins -p $period -nopsearch -dm $dm -nodmsearch -phs $phase -o prepfold_$TIMESTAMP_tube_on -noxwin ./$pulcalfil $mask_cal_on

####################################
# EDITING ON.PFD:
####################################


####################################
# REDUCTION STARTS, OFF HYDRA:
####################################
cd ..
cd reduc_off_$TIMESTAMP

# Making mask:
rfifind -time $timemask -zerodm -o mask_pul_cal_of ../$pulcalfil_of
# Folding:
prepfold -nsub $nsub -n $nbins -p $period -nopsearch -dm $dm -nodmsearch -phs $phase -o prepfold_$TIMESTAMP_tube_off -noxwin ./$pulcalfil $mask_cal_off
####################################
# EDITING OFF.PFD:
####################################

