using PropertyFunctions
using JSON
using Logging
using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using PropDicts
using FilePathsBase
using DataStructures
using PropDicts
using Tables
using TypedTables
using Optim
using FileIO
import JLD2
import HDF5

function get_settings(config)

    settings=Dict()
    settings[:energy_scale_fixed]=config["nuisance"]["energy_scale"]["fixed"]
    settings[:energy_scale_correlated]=config["nuisance"]["energy_scale"]["correlated"]
    settings[:eff_fixed]=config["nuisance"]["efficiency"]["fixed"]
    settings[:eff_correlated]=config["nuisance"]["efficiency"]["correlated"]
    settings[:bkg_only]=config["bkg_only"]

    return settings
end

function get_partitions_new(part_path::String)
    """
    Get the partition info from a JSON file and save to a Table

    """
        part_data_json = JSON.parsefile(part_path,dicttype=DataStructures.OrderedDict)

        fit_groups = part_data_json["fit_groups"]

        list_groups=collect(keys(fit_groups))
        k = keys(part_data_json["partitions"][list_groups[1]][1])
        arrays=Dict()
        for key in k
            arrays[key]=[]

        end
        arrays["fit_group"]=[]
        arrays["signal_par_name"]=[]
        arrays["bkg_par_name"]=[]
        arrays["eff_par_name"]=[]
        arrays["energy_reso_name"]=[]
        arrays["energy_bias_name"]=[]

        fit_ranges=OrderedDict()
        for fit_group in keys(part_data_json["partitions"])
            
            fit_ranges[fit_group]=part_data_json["fit_groups"][fit_group]["range"]
            for part in part_data_json["partitions"][fit_group]
                for key in k
                    append!(arrays[key],[part[key]])
                    end
                append!(arrays["fit_group"],[fit_group])
                append!(arrays["bkg_par_name"],[Symbol(part_data_json["fit_groups"][fit_group]["bkg_name"])])
            
                ## defaults to 'gaussian'
                if haskey(part_data_json["fit_groups"][fit_group], "signal_name")
                    append!(arrays["signal_par_name"], [Symbol(part_data_json["fit_groups"][fit_group]["signal_name"])])
                else
                    append!(arrays["signal_par_name"], [Symbol("gaussian")])
                end
            
                ## defaults to 'all'
                if (haskey("efficiency_group_name",part_data_json["fit_groups"][fit_group]))
                    append!(arrays["eff_par_name"],["αe_"*Symbol(part_data_json["fit_groups"][fit_group]["efficiency_group_name"])])
                else
                    append!(arrays["eff_par_name"],[:αe_all])
                end

                ## defaults to 'all'
                if (haskey("energy_scale_group_name",part_data_json["fit_groups"][fit_group]))
                    append!(arrays["energy_reso_name"],[Symbol("αr_"*part_data_json["fit_groups"][fit_group]["energy_scale_group_name"])])
                    append!(arrays["energy_bias_name"],[Symbol("αb_"*part_data_json["fit_groups"][fit_group]["energy_scale_group_name"])])

                else
                    append!(arrays["energy_reso_name"],[:αr_all])
                    append!(arrays["energy_bias_name"],[:αb_all])

                end
        
            end

        end
    
        #TODO: find a way to make this not hardcoded
        tab = Table(experiment=Array(arrays["experiment"]),
                    fit_group=Array(arrays["fit_group"]),
                    bkg_name = Array(arrays["bkg_par_name"]),
                    signal_name = Array(arrays["signal_par_name"]),
                    energy_reso_name = Array(arrays["energy_reso_name"]),
                    energy_bias_name = Array(arrays["energy_bias_name"]),
                    eff_name = Array(arrays["eff_par_name"]),
                    detector=Array(arrays["detector"]),
                    part_name=Array(arrays["part_name"]),
                    start_ts=Array(arrays["start_ts"]),
                    end_ts=Array(arrays["end_ts"]),
                    eff_tot=Array(arrays["eff_tot"]),
                    eff_tot_sigma=Array(arrays["eff_tot_sigma"]),
                    fwhm=Array(arrays["fwhm"]),
                    fwhm_sigma=Array(arrays["fwhm_sigma"]),
                    exposure=Array(arrays["exposure"]),
                    bias =Array(arrays["bias"]),
                    bias_sigma =Array(arrays["bias_sigma"]))
        return tab,fit_groups,fit_ranges
end

function get_partition_event_index(events::Array{Vector{Float64}},partitions::TypedTables.Table)::Vector{Int}
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

function get_events(event_path,partitions)::Array{Vector{Float64}}
    """
        Get the event info from a jason file and save  to a Table
    """
        @info event_path
        event_json = JSON.parsefile(event_path,dicttype=DataStructures.OrderedDict)
        events=Array{Vector{Float64}}(undef,length(partitions))
        for (idx,part) in enumerate(partitions)
            events[idx]=Vector{Float64}[]
        end
       
        for event in event_json["events"]
            found=false
            
            for (idx,part) in enumerate(partitions)
                if (part.experiment ==event["experiment"] && part.detector==event["detector"] && 
                    event["timestamp"]<part.end_ts+1 && event["timestamp"]>part.start_ts-1)
                    append!(events[idx],Vector{Float64}([Float64(event["energy"])]))
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

## sampling 
function inverse_uniform_cdf(p, fit_range)
    range_l = [arr[1] for arr in fit_range] 
    range_h = [arr[2] for arr in fit_range] 
    delta = sum(range_h .- range_l)
    a = range_l[1]
    b = range_h[end]
    
    cumulative_prob = 0.0
    cumulative_range = 0.0

    for j in 1:length(fit_range)
        interval_width = range_h[j] - range_l[j]
        interval_prob = interval_width / delta

        # if p is within the current interval range
        if cumulative_prob + interval_prob >= p
            # scale the probability to the current interval
            interval_p = (p - cumulative_prob) / interval_prob
            res = range_l[j] + interval_p * interval_width
            break
        end

        # accumulate the probability and range covered so far
        cumulative_prob += interval_prob
        cumulative_range += interval_width
    end

    # handle the p=1 case
    if p == 1
        res = range_h[end]
    end
    
    return res
end


function generate_disjoint_uniform_samples(n, fit_range)
    rands=[]
    for i in 1:n
        append!(rands,rand())
    end
    res = [inverse_uniform_cdf(rand,fit_range) for rand in rands]
    return res
end


function save_generated_samples(samples,output)
"""
Function which saves sampling results
"""
    FileIO.save(joinpath(output,"mcmc_files/samples.jld2"), Dict("samples" => samples))
    bat_write(joinpath(output,"mcmc_files/samples.h5"), samples)
end


function save_results_into_json(samples,posterior,config,output;par_names=nothing,toy_idx=nothing)
"""
Function which saves results from the fit and copies the input config (for any future need)
"""
    
    global_modes = BAT.mode(samples) 
    # a more refined estimate would be the following 
    findmode_result = bat_findmode(
        posterior,
        OptimAlg(optalg = Optim.NelderMead(), init = ExplicitInit([global_modes]))
    )
    refined_global_modes = findmode_result.result
    
    # save partitions info for nuisance parameters
    first_sample = samples.v[1]
    pars = keys(first_sample)
    nuisance_dict = Dict{String, Vector{Dict{String, String}}}()
    for par in pars
        par_entry = first_sample[par]
        
        # we do not save entries for global parameters with length==1
        if length(par_entry) != 1
            # initialize an empty array for each main key (eff, bias, res)
            nuisance_dict[string("$(par)")] = []  
            
            for idx in 1:length(par_entry) 
                xname = string("$(par)")
                if (par_names !=nothing)
                    xname = par_names[par][idx]
                end
                pattern = r"\w+\s+part\d{4}\s\w+"
                result = match(pattern, xname)
                if result !== nothing
                    tot_str = result.match
                    parts = split(tot_str)
                    inner_dict = Dict(
                        "experiment" => parts[1], 
                        "partition" => parts[2],
                        "detector" => parts[3] 
                    )
                    push!(nuisance_dict[string("$(par)")], inner_dict)
                end
            end
        end
    end
    
    ltmp = NullLogger()
    marginalized_modes=0
    with_logger(ltmp) do
        marginalized_modes = BAT.bat_marginalmode(samples).result
       end
    
    """
    unshaped_samples, f_flatten = bat_transform(Vector, samples)
    @info "Unshaped samples:", bat_report(unshaped_samples)
    
    marginalized_modes = Dict()
    parameter_samples = [s["v_1"] for s in unshaped_samples]  # Adjust for vectors if needed

    hist = fit(Histogram, parameter_samples, nbins=50)
    max_bin_idx = argmax(hist.weights)
    mode_value = hist.edges[1][max_bin_idx]

    marginalized_modes[string("v_1")] = mode_value
    println("Marginalized modes: $marginalized_modes")
    """

    mean = BAT.mean(samples)
    stddev = BAT.std(samples)
    
    ci_68 = BAT.smallest_credible_intervals(samples, nsigma_equivalent=1)
    ci_90 = BAT.smallest_credible_intervals(samples, nsigma_equivalent=1.64)
    ci_95 = BAT.smallest_credible_intervals(samples, nsigma_equivalent=2)
    ci_99 = BAT.smallest_credible_intervals(samples, nsigma_equivalent=3)
    
    quantile90 = Statistics.quantile(samples, 0.9)
    
    data = Dict(
        "mean" => mean,
        "stddev" => stddev,
        "global_modes" => global_modes,
        "refined_global_modes" => refined_global_modes,
        "marginalized_modes" => marginalized_modes,
        "ci_68" => ci_68,
        "ci_90" => ci_90,
        "ci_95" => ci_95,
        "ci_99" => ci_99,
        "quantile90" => quantile90,
        "config" => config, 
        "nuisance_partitions" => nuisance_dict
    )

    json_string = JSON.json(data,4)
    
    if toy_idx == nothing
        open(joinpath(output,"mcmc_files/fit_results.json"), "w") do file
            write(file, json_string)
        end
    else
        open(joinpath(output,"mcmc_files/fit_results_$(toy_idx).json"), "w") do file
            write(file, json_string)
        end
    end
end

function save_outputs(partitions, events, part_event_index, samples, posterior, config, output_path, fit_ranges;priors=nothing,par_names=nothing,toy_idx=nothing)
"""
Function to plot and save results, as well as inputs
"""
    if (haskey(config["bkg"],"correlated")) & (config["bkg"]["correlated"]["mode"]!="none")
        hier=true
    else
        hier=false
    end
    if (config["signal"]["prior"]=="sqrt")
        sqrt_prior=true
        s_max=config["signal"]["upper_bound"]
    else
        sqrt_prior=false
        s_max=nothing
    end
    first_sample = samples.v[1]
    free_pars = keys(first_sample) # in format (:B, :S, ...) 
    @info "... these are the parameters that were included: ", free_pars
    
    @info "... now we save samples (untouched if we do not want to overwrite)"
    if config["light_output"]==false
        if config["overwrite"]==true || !isfile(joinpath(config["output_path"],"mcmc_files/samples.h5"))
            save_generated_samples(samples, output_path)
        @info "...done!"
        end
    end
    
    @info "... now we save other useful results + config entries"
    save_results_into_json(samples, posterior, config, output_path,par_names=par_names,toy_idx=toy_idx)
    @info "...done!"

    if config["light_output"]==false
        plot_correlation_matrix(samples,output_path,par_names=par_names,toy_idx=toy_idx)
    end

    @info "... now we plot marginalized posteriors (and priors)"
    plot_marginal_distr(partitions, samples, free_pars, output_path,priors=priors,par_names=par_names,plot_config=config["plot"],s_max=s_max,sqrt_prior=sqrt_prior,hier=hier,toy_idx=toy_idx)

    if config["light_output"]==false
        @info "plot 2D posterior"
        plot_two_dim_posteriors(samples,free_pars,output_path,par_names=par_names,toy_idx=toy_idx)
        @info "...done!"
    end
    
    if config["plot"]["bandfit_and_data"] || config["plot"]["fit_and_data"]
        @info "... now we plot fit & data"
        plot_fit_and_data(partitions, events, part_event_index, samples, free_pars, output_path, config, fit_ranges, toy_idx=toy_idx)
        @info "...done!"
    end
    
end