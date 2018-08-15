#/bin/bash

#####################################
# 0) Some error checking:	
#####################################

if [ ! -f ./*$1*.fil ]; then
	echo 'ERROR: directory must contain the pulsar observation.'
	echo 'e.g.: ds256_B1240-64_020171023_104155.fil'
	exit
fi

#####################################
# 1) Pulsar information:	
#####################################

# Read pulsar name:
pulsarname=$1
pulsarname_noBJ=$(echo $pulsarname | tr -d 'B'| tr -d 'J')

# Create reduction folder:
TIMESTAMP=$(date +"%y%m%d_%H%M%S")
mkdir reduc_$TIMESTAMP

# Reading pulsar observation file name:
pulsarfil=$(ls *$pulsarname*.fil)

# Reading MJD with readfile:
MJD=$(readfile $pulsarfil | grep 'MJD' | awk '{print $5}')

# Reading observation time:
TIMESTAMP_PUL=$(readfile $pulsarfil | grep 'Obs' | awk '{print $5}')

# Computing and reading the period with period.py:
period=$(/home/observacion/data/period.py -p $pulsarname_noBJ -m $MJD | grep 'Period' | awk '{print $2}')
period=$(/home/observacion/data/period_T2.py -p $pulsarname_noBJ -m $MJD | grep 'Period' | awk '{print $2}')
period=$(bc -l <<< "$period / 1000.")

# Reading the dispersion measure from .iar file:
dm=$(less /home/observacion/iar_files/$pulsarname.iar | grep 'DM' | awk '{print $2}' | tr -d 'DM,')

#####################################
# 2) Folding parameters:	
#####################################

# Observation subchannels:
read -p 'Number of subchannels [int]: ' -e -i '32'  nsub

# Folding bins:
read -p 'Number of bins [int]: ' -e -i '512' nbins

# Pulsar period:
read -p 'Use timing? [y/n] ' -e -i 'y' timing_ans
if [ $timing_ans == "n" ] || [ $timing_ans == "n" ]; then
	echo "Then, searching period..."
	read -p 'Number of p searches [int]: ' -e -i '128' npart
	read -p 'Step for p searches [int]: ' -e -i '1' pstep
fi

# DM search:
read -p 'Dispersion measure search? [y/n] ' -e -i 'n' dmsearch
dmsearch="$(echo ${dmsearch} | tr 'A-Z' 'a-z')"
if [ $dmsearch == "y" ] || [ $dmsearch == "yes" ]; then
	dmsearch=""
	else if [ $dmsearch == "n" ] || [ $dmsearch == "no" ]; then
		dmsearch="-nodmsearch"
	else
		echo "ERROR: Wrong DM answer"
		exit
	fi
fi

# Move phase:
read -p 'Pulsar prepfold phase: ' -e -i '0.0' phs_pul

# RFI mask:
read -p 'Do pulsar mask? [y/n] ' -e -i 'n' mask_pul_ans
mask_pul_ans="$(echo ${mask_pul_ans} | tr 'A-Z' 'a-z')"
if [ $mask_pul_ans == "y" ] || [ $mask_pul_ans == "yes" ]; then
	read -p 'Time for pulsar mask integration: ' -e -i '2.0' timemask_pul
	mask_pul="-mask mask_${pulsarname}_${timemask_pul}_rfifind.mask"
	else if [ $mask_pul_ans == "n" ] || [ $mask_pul_ans == "no" ]; then
		read -p 'Reduce time interval? [y/n] ' -e -i 'n' timeint_ans
		timeint_ans="$(echo ${timeint_ans} | tr 'A-Z' 'a-z')"
		if [ $timeint_ans == "y" ] || [ $timeint_ans == "yes" ]; then
			read -p 'Start time [0:1] ' -e -i '0.0' timeint_start
			read -p 'End time [0:1] ' -e -i '1.0' timeint_end
			mask_pul="-start $timeint_start -end $timeint_end"
			else if [ $timeint_ans == "n" ] || [ $timeint_ans == "no" ]; then
			mask_pul=""
			else
				echo "ERROR: Wrong time interval answer"
				exit
			fi
		fi
	else
		echo "ERROR: Wrong pul mask answer"
	fi
fi


#####################################
# 3) REDUCTION STARTS:	
#####################################

cd reduc_$TIMESTAMP

# Making masks, if requested: 
if [ $mask_pul_ans == "y" ] || [ $mask_pul_ans == "yes" ]; then
	rfifind -time $timemask_pul -zerodm -o mask_${pulsarname}_$timemask_pul ../$pulsarfil
fi

if [ $timing_ans == "n" ] || [ $timing_ans == "n" ]; then
	prepfold -nsub $nsub -n $nbins -p $period -pstep $pstep -npart $npart -dm $dm $dmsearch -phs $phs_pul -o prepfold_$TIMESTAMP_PUL -noxwin ../$pulsarfil $mask_pul
	else
	prepfold -nsub $nsub -n $nbins -phs $phs_pul -o prepfold_$TIMESTAMP_PUL -noxwin ../$pulsarfil $mask_pul -timing /opt/pulsar_software/tempo/tzpar/$pulsarname_noBJ.par
fi

read -p 'Display prepfold output? [y/n]: ' -e -i 'y' display_pulsar
	
if [ $display_pulsar == "y" ] || [ $display_pulsar == "yes" ]; then
gv prepfold_$TIMESTAMP_PUL*.ps >/dev/null 2>&1 &
fi
cd ..

echo '########################'
echo '						  '
echo '        #         #     '
echo '        #         #     '
echo '        #         #     '
echo '        #         #     '
echo '       # #       # #    '
echo '  ### #   ##### #   ####'
echo '     #         #        '   
echo '						  '   
echo '	 puma_presto.sh: DONE '
echo '########################'