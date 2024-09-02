### main.jl
### -> gets a config.json file in input for running the Bayesian unbinned fit
###

using Pkg
Pkg.activate(".") # Activate the environment
using ArgParse

using JSON
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
    println("Reading configuration from: ", config_path)
    config = read_config(config_path)
    
    ### TO DO: if we save results somewhere (we can start from output/), it would be nice to give a name to saved files and/or stored them at specific paths - we can add this to the config file
    
    # Call the analysis function from ZeroNuFit
    ZeroNuFit.run_analysis(config)
end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
