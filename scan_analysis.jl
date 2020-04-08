using Glob, DataFrames, CSV, Plots, StatsBase, LsqFit
@. gauss(x,p) = p[1]*exp(-((x-p[2])^2)/(2*p[3]^2))

pyplot()
folders = glob("*", "output/Scan Results/")
light_yields = []
abs_vals = []
abs_val_input = 10:5:25
light_yield_input = 3400:200:6800

peaks = zeros(Float64, length(abs_val_input), length(light_yield_input))
#sigma of 32
for folder in values(folders)
           #println(folder)
           files = glob("*_t*.csv", folder)
           abs = parse(Int,files[1][64:65])
           ly = parse(Int, files[1][72:75])
           #println(abs ,"  ", ly)
           df = DataFrame(N_Trigger = Int64[], EDep_Trigger = Float32[], EDep_PEN = Float32[], N_Left = Int64[], N_Right = Int64[], N_Bottom = Int64[], N_Front = Int64[], N_Back = Int64[])
           for file in values(files)
               data = CSV.read(file, comment="#", header = ["N_Trigger","EDep_Trigger", "EDep_PEN", "N_Left", "N_Right", "N_Bottom", "N_Front", "N_Back"])
               append!(df, data)
           end
           df.N_total = df.N_Left .+ df.N_Right  .+ df.N_Front .+ df.N_Back #.+ df.N_Bottom

           h_nphotons = fit(Histogram, df.N_total, 0:20:6000)
           peak_bin = findmax(h_nphotons.weights[10:end])[2]+9
           #println(peak_bin)
           gauss_fit = curve_fit(gauss, midpoints(h_nphotons.edges[1])[peak_bin-10:peak_bin+10], h_nphotons.weights[peak_bin-10:peak_bin+10], [maximum(h_nphotons.weights), midpoints(h_nphotons.edges[1])[peak_bin],10])

           # plot(h_nphotons, st=:step, xlabel = "N Photons", ylabel = "Counts / 20 Photons");
           # plot!(midpoints(h_nphotons.edges[1])[peak_bin-10:peak_bin+10],x->gauss(x,gauss_fit.param), label = "$(gauss_fit.param[2])");
           # savefig("$(folder)/output_plot.pdf");

           push!(abs_vals, abs)
           push!(light_yields, ly)
           abs_bin = findfirst(x->x==abs, abs_val_input)

           ly_bin =  findfirst(x->x==ly, light_yield_input)
           #println(abs_bin ,"  ", ly_bin)
           peaks[abs_bin, ly_bin] = gauss_fit.param[2]
end
heatmap(light_yield_input, abs_val_input, peaks, xlabel = "Input Light Yield [Photons / MeV]", ylabel = "Input Absortion length @ 450 nm [mm]", size = (400,400))
savefig("heatmap_test.pdf")
