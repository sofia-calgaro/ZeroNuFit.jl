using LegendDataManagement
using PropertyFunctions
using JSON
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


function save_results(samples,output)
    FileIO.save(joinpath(output,"mcmc_files/results.jld2"), Dict("samples" => samples))
    bat_write(joinpath(output,"mcmc_files/results.h5"), samples)
end


