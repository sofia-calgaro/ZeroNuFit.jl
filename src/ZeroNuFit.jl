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
    
    @info"Let's retrieve some partitions ..."
    partitions = []
    for part_path  in config["partitions"]
        append!(partitions,[get_partitions_new(part_path)])
    end
    @info "... load events"
    events = []
    for (event_path,part) in zip(config["events"],partitions)
        append!(events,[get_events(event_path,part)])
    end
    @debug "... extracted events:", events

    @info "get which partitions have events"
    part_event_index = get_partition_event_index(events[1],partitions[1])
    
    # check if you want to overwrite the fit; if no results are present, then fit data
    if config["overwrite"] == true || !isfile(joinpath(config["output_path"],"mcmc_files/samples.h5"))
        @info "... now we run a fit"
        if config["overwrite"] == true
            @info "OVERWRITING THE PREVIOUS FIT!"
        end
        samples = run_fit_over_partitions(partitions[1],events[1],part_event_index,config=config,stat_only=config["stat_only"]) 
    # load the already present fit
    else
        @info "... we load already existing fit results"
        samples = bat_read(joinpath(config["output_path"],"mcmc_files/samples.h5")).result
    end
    
    @info bat_report(samples)
    
    # save results
    save_outputs(partitions[1], samples, config)
    
    return 
end

end 
