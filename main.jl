include("control.jl")
#include("analysis.jl")
using CSV

abs_files = ["400nm-Fixed-20mm", "450nm-Fixed-20mm"]
abs_vals = 20:5:25
light_yields = 3400:200:6800
n_events = 1000000
#n_events = 50

input_folder = "input_files/"
top_dir = pwd()

for abs_val in abs_vals
    for ly in light_yields

        input_data = CSV.read(input_folder*"Exp4.csv")
        scale_factor =  abs_val / input_data[2][50]
        input_data[2] = input_data[2][:].*scale_factor
        CSV.write(input_folder*"450nm-Fixed-$(abs_val)mm.csv",input_data)
        cd("build/")
        script = "
        /PEN/det/setDetName 450nm-Fixed-abs_$(abs_val)mm-LY-$ly
        /PEN/det/setAbsFile 450nm-Fixed-$(abs_val)mm
        /PEN/det/setLY $ly
        /run/beamOn $n_events"

       open("test.mac","w") do f
           write(f,script)
       end
        run(`./PEN -m test.mac -t 24`)
        cd(top_dir)
    end
end
