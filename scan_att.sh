#!/bin/bash
declare -a position=(-25. -20. -15. -10. -5. 0. 5. 10. 15. 20. 25.)
#declare -a LY=(5500. 6500. 7500.)
declare -a LY=(5500.)
#declare -a attenuation=(10. 15. 20. 25.)
declare -a attenuation=(10.)
cd ./build/
for i in "${position[@]}"
do
	echo item: $i
	for j in "${LY[@]}"
	do
		echo item: $j
		for k in "${attenuation[@]}"
		do
			echo item: $k
			rm -rf Bi207_scan.mac
			echo "/PEN/det/setDetectorType 1
/PEN/det/setCollimatorPositionX $i mm
/PEN/det/setLY $j
/PEN/det/setABS $k
/PEN/gun/sourcePositionX $i mm
/PEN/gun/sourceType 1

/run/beamOn 10000000" >> Bi207_scan.mac
			./PEN -m Bi207_scan.mac -t 24
			sleep 1 
		done
	done
done
















































