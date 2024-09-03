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
    
    @info "... and now we run a fit"
    samples_uniform = run_fit_over_partitions(partitions[1],events[1],func=fit_function_linear,config=config,stat_only=config["stat_only"]) 
    println(samples_uniform)
    println(bat_report(samples_uniform))
    
    @info "and we plot results"
    energies = []
    for (idx_k, part_k) in enumerate(partitions[1])
        if events[1][idx_k] != Any[]
            for energy in events[1][idx_k]
                append!(energies, events[1][idx_k])
            end
        end
    end
    hist_data = append!(Histogram(1930:2:2190), energies)
    make_plots(samples_uniform, (:B, :S), hist_data, fit_function_linear, config["output_path"])
    
    @info "and we save fit results"
    save_results(samples_uniform,config["output_path"])
    @info "...done!"
    
    return 
end

end 
