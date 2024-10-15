# src/ZeroNuFit.jl
module ZeroNuFit

include("fitting.jl")
include("plotting.jl")
include("utils.jl")
using JSON

export run_analysis
export retrieve_real_fit_results



function get_partitions_events(config::Dict{String, Any})
    @info "You entered into src/ZeroNuFit.jl"
    
    @info"Let's retrieve some partitions ..."
    partitions = nothing
    first=true
    fit_ranges=nothing
    for part_path  in config["partitions"]

        part_tmp,fit_groups,fit_range =get_partitions_new(part_path) 
        if (first)
            partitions=part_tmp
            first=false
            fit_ranges=fit_range
        else
            partitions=vcat(partitions,part_tmp)
            merge!(fit_ranges,fit_range)
        end
    end
    display(partitions)
    @info "... load events"
    events_multi = []
    for event_path in config["events"]
        append!(events_multi,[get_events(event_path,partitions)])
    end

    events=Array{Vector{Float64}}(undef,length(partitions))
    for i in 1:length(partitions)
        
        arr_tmp =Vector{Float64}()
        for sub in events_multi
            if (sub[i]!=Float64[])
                
                append!(arr_tmp,sub[i])
                end
        end
         
        events[i]=arr_tmp
    end
    @debug events

    @info "get which partitions have events"
    part_event_index = get_partition_event_index(events,partitions)

    return part_event_index,events,partitions,fit_ranges

end
# function to run the unbinned fit
function run_analysis(config::Dict{String, Any};output_path::String, toy_idx=nothing)
"""
Function which handeles running analysis
Parameters:
----------
    config::Dict{String,Any} the fit configuration
    output_path::String (keyword) the path to the output files folder
"""

   
    part_event_index,events,partitions,fit_ranges= get_partitions_events(config)
    # check if you want to overwrite the fit; if no results are present, then fit data
    if config["overwrite"] == true || !isfile(joinpath(config["output_path"],"mcmc_files/samples.h5"))
        @info "... now we run a fit"

        if config["overwrite"] == true
            @info "OVERWRITING THE PREVIOUS FIT!"
        end

        samples,prior,par_names = run_fit_over_partitions(partitions,events,part_event_index,config,fit_ranges) 
        @info "fit ran succesfully"
    else
        @info "... we load already existing fit results"
        samples = bat_read(joinpath(config["output_path"],"mcmc_files/samples.h5")).result
        prior,_,_,par_names = get_stat_blocks(partitions,events,part_event_index,fit_ranges,config=config,bkg_only=config["bkg_only"]) 
    end

    # let's save
    @info samples
    @info bat_report(samples)
    _,_,posterior,_ = get_stat_blocks(partitions,events,part_event_index,fit_ranges,config=config,bkg_only=config["bkg_only"]) 
    save_outputs(partitions, events, part_event_index, samples, posterior, config, output_path, fit_ranges, priors=prior,par_names=par_names, toy_idx=toy_idx)
    
    return 
    
end


function retrieve_real_fit_results(config::Dict{String, Any})
"""
Function which handeles generating of fake data
Parameters:
----------
    config::Dict{String,Any} the fit configuration
    output_path::String (keyword) the path to the output files folder
"""

    @info "You entered into src/ZeroNuFit.jl"
    
    @info"Let's retrieve some partitions ..."
    partitions = nothing
    first=true
    @info config["partitions"]
    for part_path  in config["partitions"]

        part_tmp,fit_groups =get_partitions_new(part_path) 
        if (first)
            partitions=part_tmp
            first=false
        else
            partitions=vcat(partitions,part_tmp)
        end
    end
    display(partitions)
    @info "... load events"
    events_multi = []
    for event_path in config["events"]
        append!(events_multi,[get_events(event_path,partitions)])
    end

    events=Array{Vector{Float64}}(undef,length(partitions))
    for i in 1:length(partitions)
        
        arr_tmp =Vector{Float64}()
        for sub in events_multi
            if (sub[i]!=Float64[])
                
                append!(arr_tmp,sub[i])
                end
        end
         
        events[i]=arr_tmp
    end
    @debug events

    @info "get which partitions have events"
    part_event_index = get_partition_event_index(events,partitions)

    # let's retrieve the old data
    samples = bat_read(joinpath(config["output_path"],"mcmc_files/samples.h5")).result

    return samples, partitions, part_event_index

end 

end