### main.jl
### -> gets a config.json file in input for running the Bayesian unbinned fit
###

using Logging
using Pkg
Pkg.activate(".") # Activate the environment
using ArgParse

using JSON
using FilePathsBase
# load the script to run the analysis
include("src/ZeroNuFit.jl")
using .ZeroNuFit


# read JSON configuration file
function read_config(file_path::String)
    json_string = read(file_path, String)
    config = JSON.parse(json_string)
    return config
end

# process parsed arguments for the main function
function get_argparse()
    settings = ArgParseSettings(prog="LEGEND ovbb Bayesian unbinned fit",
                            description="",
                            commands_are_required = true)
    @add_arg_table settings begin
        "--config", "-c"
            help = "path to config file"
            arg_type = String
            required = true
    end
    
    parse_args(settings)
    return parse_args(settings)
end

function main()
    
    # read parsed arguments
    parsed_args = get_argparse()
    # read config path
    config_path = parsed_args["config"]
    @info "Reading configuration from: $config_path"
    config = read_config(config_path)
    
    # load the output path and create the neccesary
    output_path = config["output_path"]

    for dir in ["$output_path/","$output_path/plots/","$output_path/mcmc_files/","$output_path/logs/"]
        if !isdir(dir)
            mkpath(dir)
        end
    end
   

    # Call the analysis function from ZeroNuFit
    ZeroNuFit.run_analysis(config,output_path=output_path)

end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
