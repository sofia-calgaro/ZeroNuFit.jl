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
    

    
    @info"Let's try retrieving some partitions ..."
    partitions = []
    for part_path  in config["partitions"]
        append!(partitions,[get_partitions_new(part_path)])
    end
    @info "... load events"
    events = []
    for (event_path,part) in zip(config["events"],partitions)
        append!(events,[get_events(event_path,part)])
    end
    @debug "extracted events", events
    @info "...done!"
    
    @info "and now we run a fit"
    samples_uniform = run_fit_over_partitions(partitions[1],events[1],func=fit_function_uniform,config=config,stat_only=config["stat_only"]) 
    #bat_report(samples_l200_uniform)

    @info "...done!"
    
    return 
end

end 
