# src/ZeroNuFit.jl
module ZeroNuFit

include("fitting.jl")
include("plotting.jl")
include("utils.jl")
using JSON

export run_analysis

# function to run the unbinned fit
function run_analysis(config::Dict{String, Any})
    println("You entered into sr/ZeroNuFit.jl")
    println("Starting the analysis...")
    println("- running analysis with the following configuration:")
    println(config)
    
    l200 = [1953.1427, 1955.2213, 1974.731, 1996.4917, 2016.76, 2040.262, 2095.7217]
    
    #samples_l200_uniform = run_fit(l200,fit_function_uniform,distprod(b=0..20.))
    #bat_report(samples_l200_uniform)
    get_partitions(config)
    
    println("...done!")
    return # fit_results
end

end 
