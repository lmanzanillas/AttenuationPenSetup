"""
    run_sim()
    Start PEN simulation run of given macro file with given number of cores.
    Default to run a single threaded test.mac
"""

function run_sim(;macro_file::String="test.mac", ncores::Int64=1)
    cd("build/")
    run(`./PEN -m $macro_file -t $ncores`)
    cd("..")
end
