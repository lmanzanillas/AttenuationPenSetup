
abs_vals = 10.:5.:25.
collimator_position = 5.:5.:25.
light_yields = 4500.:1000.:7500.
n_events = 10000000

top_dir = pwd()

for abs_val in abs_vals
    for ly in light_yields
	for pos in collimator_position
       	 	cd("build/")
        	script ="
	        /PEN/det/setDetectorType 1
		/PEN/det/setCollimatorPositionX $pos mm
        	/PEN/det/setLY $ly
		/PEN/det/setABS $abs_val
		/PEN/det/setSigAlpha 0.7
		/PEN/gun/sourcePositionX $pos mm
		/PEN/gun/sourceType 1

        	/run/beamOn $n_events"

	        open("Bi_scan.mac","w") do f
        	   write(f,script)
       		end
        	run(`./PEN -m Bi_scan.mac -t 24`)
        	cd(top_dir)
	end
    end
end
