using CSV, DataFrames, StatsBase, Glob, Plots, LsqFit

@. model(x, p) = (1 / sqrt(2*pi*p[2]^2)) * exp(-x-p[2]^2/(2*p[2]^2))

Plots.pyplot()

"""
    run_analysis(folder::String)
    Read and plot simulation data from folder.
    Saves plot of energy deposit and detected light output from PEN.
"""

function run_analysis(folder::String)
    files = glob("output/$folder/*_t*.csv")
    df = DataFrame(N_Trigger = Int64[], EDep_Trigger = Float32[], EDep_PEN = Float32[], N_Left = Int64[], N_Right = Int64[], N_Bottom = Int64[], N_Front = Int64[], N_Back = Int64[])
    for file in values(files)
        data = CSV.read(file, comment="#", header = ["N_Trigger","EDep_Trigger", "EDep_PEN", "N_Left", "N_Right", "N_Bottom", "N_Front", "N_Back"])
        append!(df, data)
    end
    df.N_total = df.N_Left .+ df.N_Right .+ df.N_Bottom .+ df.N_Front .+ df.N_Back

    h_nphotons = fit(Histogram, df.N_total, 0:20:6000)
    peak_bin = findmax(h_nphotons.weights[30:end])[2]+29
    photon_plt = plot(h_nphotons, st=:step, label = "PEN", legend = :topright, xlabel = "Number of Photons", ylabel = "Events / 10 photons")

    h_edep_pen = fit(Histogram, df.EDep_PEN, 0:10:1200)
    pen_plt = plot(h_edep_pen, st=:step, label = "Foil", legend = :topleft, xlabel = "Energy in PEN [keV]")

    vline!(photon_plt, [midpoints(0:20:6000)[peak_bin]], label = "$(midpoints(0:20:6000)[peak_bin])")
    combi_energy = plot(xlabel = "Energy [keV]", ylabel = "Events / 20 keV")
    plot!(combi_energy, h_edep_pen, st=:step, label = "PEN")
    combined_plt = plot(layout = (2,1), size = (800,800), combi_energy, photon_plt)
    savefig(combined_plt, "output/Plots/pen_sim_$folder.pdf")
    return photon_plt
end
