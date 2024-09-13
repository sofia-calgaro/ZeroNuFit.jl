### main.jl
### -> gets a config.json file in input for running the Bayesian unbinned fit
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
    
    # read parsed arguments
    @info "running using ",Base.Threads.nthreads()," threads"
    parsed_args = get_argparse()
    
    # read N toy index
    toy_idx = parsed_args["index_toy"]

    # read config path
    config_path = parsed_args["config"]
    @info "Reading configuration from: $config_path"
    config = read_config(config_path)
    
    # load the output path and create the neccesary
    output_path = config["path_to_fit"]
    
    config_real_data = read_config(joinpath(config["path_to_fit"], "mcmc_files/fit_results.json"))["config"]
    # we add/overwrite an option for light saving
    config_real_data["light_output"] = true

    for dir in ["$output_path/sensitivity","$output_path/sensitivity/fake_data/","$output_path/sensitivity/plots/","$output_path/sensitivity/mcmc_files/","$output_path/sensitivity/logs/"]
        if !isdir(dir)
            mkpath(dir)
        end
    end
   
    set_logger(config_real_data,"$output_path/sensitivity",toy_idx=toy_idx)
    
    # let's retrieve input for the fake generation of data (JUST ONCE!)
    samples, partitions, part_event_index = retrieve_real_fit_results(config_real_data)
    
    # now let's generate and fit data! How many times? As N_toys
    fake_data = generate_data(samples,partitions,part_event_index,best_fit=config["best_fit"],nuis_prior=config_real_data["nuisances"]["prior"],bkg_only=config["bkg_only"],seed=config["seed"])

    # define a new path for the events (where we will save everything)
    config_real_data["events"] = ["$output_path/sensitivity/fake_data/fake_data$toy_idx.json"]
    # save fake data there
    open(config_real_data["events"][1], "w") do file
        JSON3.write(file, fake_data)
    end

    # fit again (and plot fake data + fit result)
    config_real_data["plot"]["fit_and_data"] = true
    ZeroNuFit.run_analysis(config_real_data,output_path="$output_path/sensitivity",toy_idx=toy_idx)
        
end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
