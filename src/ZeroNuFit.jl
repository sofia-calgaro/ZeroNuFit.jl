# src/ZeroNuFit.jl
module ZeroNuFit

using JSON

export run_analysis

function run_analysis(config::Dict{String, Any})
    println("You entered into sr/ZeroNuFit.jl - starting the analysis ...")
    println("Running analysis with the following configuration:")
    println(config)
    
    # ... here the analysis will be implemented ...
    
    println("...done!")
end

end 
