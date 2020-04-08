include("control.jl")
include("analysis.jl")

abs_files = ["400nm-Fixed-20mm", "450nm-Fixed-20mm"]
light_yields = [3400, 6800]
n_events = 5000000

cd("build/")
for abs_file in abs_files
    for ly in light_yields
    script = "/PEN/det/setDetName $(abs_file[1:end-4])-LY-$ly
       /PEN/det/setAbsFile $abs_file
       /PEN/det/setLY $ly
       /run/beamOn $n_events"

       open("test.mac","w") do f
           write(f,script)
       end
        run(`./PEN -m test.mac -t 24`)
    end
end
include("analysis.jl")
fixed_450_6800 = run_analysis("Bi207-20-03-2020 07-04-05")
fixed_450_3400 = run_analysis("Bi207-20-03-2020 04-02-06")
fixed_400_6800 = run_analysis("Bi207-19-03-2020 23-47-45")
fixed_400_3400 = run_analysis("Bi207-19-03-2020 20-10-46")
a = plot(fixed_450_6800, title = "450nm Fixed, 6800 photon/MeV")
b = plot(fixed_450_3400, title = "450nm Fixed, 3400 photon/MeV")
c = plot(fixed_400_6800, title = "400nm Fixed, 6800 photon/MeV")
d = plot(fixed_400_3400, title = "400nm Fixed, 3400 photon/MeV")

merge = plot(a, b, c, d, layout = (2,2), size = (800, 800))
savefig(merge, "sim_output.pdf")
