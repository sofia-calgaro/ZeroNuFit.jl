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

function get_partitions(config::Dict{String, Any}, printflag=false)
    # retrieve metadata path from the config; if nothing, load the online metadata
    metapath = config["meta_path"]
    if metapath != nothing
        meta = LegendDataManagement.AnyProps(metapath)
    else
        meta = LegendDataManagement.AnyProps() ### TO DO: test if it works
    end

    partitions = meta.datasets.ovbb_partitions_pars
    function printdb(db)
        if db isa PropDict
            for (key,value) in db
                println(key)
                printdb(value)
            end
        else
            println("\t",db)
            println()
        end
    end
    
    # loop over detectors
    for (detector, detdata) in partitions
        if detector == :default
            continue
        end

        if printflag == true
            println("detector = $detector")
            println(">>>>>>>>>>> ")
        end
        
        # apply defaults
        if partitions.default isa PropDicts.MissingProperty 
            detdata_merge = copy(detdata)
        else
            new = partitions.default
            detdata_merge = merge(new,copy(detdata))
        end
        
        # loop over partitions for this detector
        for (partition, pardata) in detdata_merge
            if partition == :default
                continue
            end

            if detdata_merge.default isa PropDicts.MissingProperty 
                pardata_merge = copy(pardata)
            else
                new = detdata_merge.default
                pardata_merge = merge(new,copy(pardata))
            end

            if printflag == true
                println("for partition $partition ")
                printdb(pardata_merge)
                println("\n")
            end

        end
    end
end