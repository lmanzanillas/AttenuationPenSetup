#!/bin/bash
#declare -a position=(-25. -20. -15. -10. -5. 0. 5. 10. 15. 20. 25.)
#declare -a position=(-10. -5. 0. 5. 10. 15. 20. 25.)
#declare -a LY=(5500. 6500. 7500.)
#declare -a LY=(5500.)
#declare -a attenuation=(10. 15. 20. 25.)
declare -a position=(0.)
declare -a SigmaAlpha=(0.5 0.9)
declare -a Reflectivity=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
cd ./build/
for i in "${position[@]}"
do
	echo item: $i
	for j in "${SigmaAlpha[@]}"
	do
		echo item: $j
		for k in "${Reflectivity[@]}"
		do
			echo item: $k
			rm -rf Bi207_scan.mac
			echo "#Scan with Bi
/PEN/det/setDetectorType 2
/PEN/det/setNTargetSamples 1
/PEN/det/setTargetThickness 10. mm
/PEN/det/setCollimatorPositionX $i mm
/PEN/det/setTargetMat PVT_structure
/PEN/det/setSigAlpha $j
/PEN/det/setPMTReflectivity $k
/PEN/gun/sourcePositionX $i mm
/PEN/gun/sourceType 1

/run/beamOn 1000000" >> Bi207_scan.mac
			./PEN -m Bi207_scan.mac -t 16
			sleep 1 
		done
	done
done
















































