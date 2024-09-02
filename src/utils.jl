using LegendDataManagement
using PropertyFunctions
using JSON
using PropDicts

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