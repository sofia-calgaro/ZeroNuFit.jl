using LegendDataManagement
using PropertyFunctions
using JSON
using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using PropDicts
using FilePathsBase
using DataStructures
using PropDicts
using Tables
using TypedTables
using LegendDataManagement
using FileIO
import JLD2
import HDF5

function get_partitions_new(part_path::String)
    """
    Get the partition info from a jason file and save  to a Table

    """
        part_data_json = JSON.parsefile(part_path,dicttype=DataStructures.OrderedDict)
        k = keys(part_data_json[1])
        arrays=Dict()
        for key in k
            arrays[key]=[]
        end
        for part in part_data_json
            for key in k
                append!(arrays[key],[part[key]])
                end
    
        end

        
        #TODO: find a way to make this not hardcoded
        tab = Table(detector=Array(arrays["detector"]),
                    start_ts=Array(arrays["start_ts"]),
                    end_ts=Array(arrays["end_ts"]),
                    eff_tot=Array(arrays["eff_tot"]),
                    eff_tot_sigma=Array(arrays["eff_tot_sigma"]),
                    fwhm=Array(arrays["fwhm"]),
                    fwhm_sigma=Array(arrays["fwhm_sigma"]),
                    exposure=Array(arrays["exposure"]),
                    bias =Array(arrays["bias"]),
                    bias_sigma =Array(arrays["bias_sigma"]))
    
        return tab
end

function get_partition_event_index(events,partitions)
"""
gets an object descirbing if a partiton has an event and giving them indexs
This creates a vector where
V[i]=0 if partition i has no events
V[i]=idx if partition i has events

where the index counts the number of partitions with index<=i with, 
events and corresponds to the index of the parameters.

"""

    output = Vector{Int}(undef,length(partitions))
    counter=1
    for (idx,part) in enumerate(partitions)
        if (events[idx] != Any[])
            output[idx]=counter
            counter+=1
        else
            output[idx]=0
        end
    end
    return output
end

function get_events(event_path,partitions)
    """
        Get the event info from a jason file and save  to a Table
    """
    
        event_json = JSON.parsefile(event_path,dicttype=DataStructures.OrderedDict)
        events=[]
        for (idx,part) in enumerate(partitions)
            append!(events,[[]])
        end
       
        for event in event_json["events"]
            found=false
            
            for (idx,part) in enumerate(partitions)
                if (part.detector==event["detector"] && event["timestamp"]<part.end_ts && event["timestamp"]>part.start_ts)
                    append!(events[idx],event["energy"])
                    found=true
                end
                    
            end
            if (found==false)
                @error event "has no partition"
                exit(-1)
            end
        end
        return events
        
end


function save_generated_samples(samples,output)
"""
Function which saves sampling results
"""
    FileIO.save(joinpath(output,"mcmc_files/samples.jld2"), Dict("samples" => samples))
    bat_write(joinpath(output,"mcmc_files/samples.h5"), samples)
end


function save_results_into_json(samples,config,output)
"""
Function which saves results from the fit and copies the input config (for any future need)
"""
    
    global_modes = BAT.mode(samples) # quick estimate
    @info "Global modes: ", global_modes
    marginalized_modes = BAT.bat_marginalmode(samples).result
    @info "Marginalized modes: ", marginalized_modes
    mean = BAT.mean(samples)
    @info "Mean: ", mean
    stddev = BAT.std(samples)
    @info "Stddev: ", stddev
    
    ##a more refined estimate would be the following 
    #using Optim
    #findmode_result = bat_findmode(
    #    samples,
    #    OptimAlg(optalg = Optim.NelderMead(), init = ExplicitInit([samples_mode]))
    #)
    #fit_par_values = findmode_result.result
    
    """
    # to be fixed
    ci_68 = credible_interval(samples, 0.68)
    ci_90 = credible_interval(samples, 0.90)
    ci_95 = credible_interval(samples, 0.95)
    
    credible_intervals = Dict()
    for par in keys(ci_68)
        credible_intervals[par] = Dict(
            "68%" => ci_68[par],
            "90%" => ci_90[par],
            "95%" => ci_95[par]
        )
    end
    """

    data = Dict(
        "mean" => mean,
        "stddev" => stddev,
        "global_modes" => global_modes,
        "marginalized_modes" => marginalized_modes,
        #"credible_interval" => credible_intervals,
        "config" => config
    )

    json_string = JSON.json(data)
    
    open(joinpath(output,"mcmc_files/fit_results.json"), "w") do file
        write(file, json_string)
    end

end

function save_outputs(samples, config)
"""
Function to plot and save results, as well as inputs
"""
    
    output_path = config["output_path"]
    free_pars = (:B, :S) # paramnames(samples) why it does not work? TO DO
    @info "... these are the parameters that were included: ", free_pars
    
    @info "... now we plot results"
    make_plots(samples, free_pars, output_path)
    @info "...done!"
    
    @info "... now we save samples (untouched if we do not want to overwrite)"
    if config["overwrite"] == false
        save_generated_samples(samples, output_path)
    @info "...done!"
    end
    
    @info "... now we save other useful results + config entries"
    save_results_into_json(samples, config, output_path)
    @info "...done!"
    
end