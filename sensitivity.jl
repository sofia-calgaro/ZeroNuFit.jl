### sensitivity.jl
### -> generate or retrieve toy data and run an unbinned Bayesian fit
###

# load the script to run the analysis
using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
using JSON3
include("main.jl")
include("src/likelihood.jl")
include("src/utils.jl")


function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--config", "-c"
            help = "Path to config file"
            arg_type = String
            required = true
        "--index_toy", "-i"
            help = "Index of the toy to generate and run the fit. If nothing is provided and there is a path to already generated toy JSON data in the configuration JSON file in input, then you will retrieve those already present fake data. "
            arg_type = Int
            default = nothing
            required = true
        "--path_to_toys", "-p"
            help = "Path to fake data folder with already existing toys (not necessarily required)"
            arg_type = String
            default = nothing
            required = false
    end

    parsed_args = parse_args(s)
    path_to_toys = parsed_args["path_to_toys"]
    
    # read parsed arguments
    @info "running using ",Base.Threads.nthreads()," threads"
    
    # read N toy index
    toy_idx = parsed_args["index_toy"]

    # read config path
    config_path = parsed_args["config"]
    @info "Reading configuration from: $config_path"
    config = read_config(config_path)
    
    # load the output path and create the neccesary
    output_path = config["path_to_fit"]
    saving_folder = get(config, :saving_folder, "sensitivity") 
    
    config_real_data = read_config(joinpath(config["path_to_fit"], "mcmc_files/fit_results.json"))["config"]
    # we add/overwrite an option for light saving
    config_real_data["light_output"] = true

    for dir in ["$output_path/$saving_folder","$output_path/$saving_folder/fake_data/",
                "$output_path/$saving_folder/plots/","$output_path/$saving_folder/mcmc_files/",
                "$output_path/$saving_folder/logs/"]
        if !isdir(dir)
            mkpath(dir)
        end
    end
   
    set_logger(config_real_data,"$output_path/$saving_folder",toy_idx=toy_idx)
    
    # we generate+fit a toy spectrum
    if path_to_toys == nothing
        @info "You'll generate new toys!"
        # let's retrieve input for the fake generation of data (JUST ONCE!)
        samples, partitions, part_event_index = retrieve_real_fit_results(config_real_data)

        # get fit ranges
        fit_ranges=nothing
        first=true
        for part_path  in config_real_data["partitions"]
            _,_,fit_range =get_partitions_new(part_path) 
            if (first)
                first=false
                fit_ranges=fit_range
            else
                merge!(fit_ranges,fit_range)
            end
        end

        # now let's generate and fit data! How many times? As N_toys
        settings=get_settings(config_real_data)
        fake_data = generate_data(samples,partitions,part_event_index,settings,fit_ranges,best_fit=config["best_fit"],seed=config["seed"])

        # define a new path for the events (where we will save everything)
        config_real_data["events"] = ["$output_path/$saving_folder/fake_data/fake_data$toy_idx.json"]
        # save fake data there
        open(config_real_data["events"][1], "w") do file
            JSON3.write(file, fake_data)
        end
        
    # we retrieve+fit a toy spectrum
    else
        @info "You'll retrieve already existing toys!"
        if !ispath(path_to_toys)
            @info "Path to toy data does not exist! Exit here."
            exit()
        end
        # save fake data again (in future, we can introduce a symlink to already existing files)
        cp("$path_to_toys/fake_data$toy_idx.json", "$output_path/$saving_folder/fake_data/fake_data$toy_idx.json", force=true)
        config_real_data["events"] = ["$output_path/$saving_folder/fake_data/fake_data$toy_idx.json"]
        @info "We copied $path_to_toys/fake_data$toy_idx.json into $output_path/$saving_folder/fake_data/fake_data$toy_idx.json"
    end
    
    # enable plotting of fit over fake data
    config_real_data["plot"]["fit_and_data"] = true
    
    # include signal in the fit (and check if previously S=0 or not - we just raise a warning)
    if config_real_data["bkg_only"] == false
        @info "*** The fit over real data was done with B!=0 and S!=0 ***"
    else
        @info "*** The fit over real data was done with B!=0 and S=0 ***"
    end
    config_real_data["bkg_only"] = false
    @info "...we now set the fit option over fake data to B!=0 and S!=0"
    
    if config["low_stat"]==true
        @info "we set the MCMC statistics to 10^5 iterations and 5 chains"
        config_real_data["bat_fit"]["nsteps"] = 1e5
        config_real_data["bat_fit"]["nchains"] = 5
    end

    # fit fake data
    ZeroNuFit.run_analysis(config_real_data,output_path="$output_path/$saving_folder",toy_idx=toy_idx)
        
end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
