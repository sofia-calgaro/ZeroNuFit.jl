# src/ZeroNuFit.jl
module ZeroNuFit

include("fitting.jl")
include("plotting.jl")
include("utils.jl")
using JSON

export run_analysis

# function to run the unbinned fit
function run_analysis(config::Dict{String, Any};output_path::String)
"""
Function which handeles running analysis
Parameters:
----------
    config::Dict{String,Any} the fit configuration
    output_path::String (keyword) the path to the output files folder
"""

    @info "You entered into src/ZeroNuFit.jl"
    
    @info"Let's retrieve some partitions ..."
    partitions = nothing
    first=true
    for part_path  in config["partitions"]

        part_tmp,fit_groups =get_partitions_new(part_path) 
        if (first)
            partitions=part_tmp
            first=false
        else
            partitions=vcat(partitions,part_tmp)
        end
    end
    @info display(partitions)
    @info "... load events"
    events = []
    for event_path in config["events"]
        append!(events,get_events(event_path,partitions))
    end
    @debug "... extracted events:", events

    @info "get which partitions have events"
    part_event_index = get_partition_event_index(events,partitions)
    
    # check if you want to overwrite the fit; if no results are present, then fit data
    if config["overwrite"] == true || !isfile(joinpath(config["output_path"],"mcmc_files/samples.h5"))
        @info "... now we run a fit"
        if config["overwrite"] == true
            @info "OVERWRITING THE PREVIOUS FIT!"
        end
        samples = run_fit_over_partitions(partitions,events,part_event_index,config=config,stat_only=config["stat_only"]) 
        @info "fit ran succesfully"
    else
        @info "... we load already existing fit results"
        samples = bat_read(joinpath(config["output_path"],"mcmc_files/samples.h5")).result
    end
    
    @info bat_report(samples)
    
    save_outputs(samples, config)
    
    return 
end

end 
