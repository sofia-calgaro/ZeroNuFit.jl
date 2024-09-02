# src/ZeroNuFit.jl
module ZeroNuFit

include("fitting.jl")
include("plotting.jl")
include("utils.jl")
using JSON

export run_analysis

# function to run the unbinned fit
function run_analysis(config::Dict{String, Any};output_path::String)

    @info "You entered into src/ZeroNuFit.jl"
    @info "Starting the analysis..."
    @info "- running analysis with the following configuration:"
    
    
    
    @info "we define some legend data: $l200"
    
    @info"Let's try retrieving some partitions ..."
    partitions = []
    for part_path  in config["partitions"]
        append!(partitions,get_partitions_new(part_path))
    end
    events = []
    for (event_path,part) in zipped(config["events"],partitions)
        append!(events,get_events(event_path,part))
    end
    @info "...done!"
    
    @info "and now we run a fit:"
    
    ### TO DO: some of these specifications must go in the config file
    # this is a test function
    
    samples_l200_uniform = run_fit_over_partitions(partitions,l200,fit_function_uniform,prior) 
    #bat_report(samples_l200_uniform)

    @info "...done!"
    
    return 
end

end 
