### main.jl
### -> gets a config.json file in input for running the Bayesian unbinned fit
###


using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate() # instantiate the environment
using ArgParse
using Logging, LoggingExtras
using JSON
using FilePathsBase
# load the script to run the analysis
include("src/ZeroNuFit.jl")
using .ZeroNuFit


function set_logger(config::Dict,output_path::String;toy_idx=nothing)
"""
Function which sets the logging for the program
Parameters
----------
    config::Dict the fit config
    output_path::String path to save the logs to
"""
    if ("debug" in keys(config) && config["debug"]==true)
        terminal_log=global_logger(ConsoleLogger(stderr, LogLevel(Debug)))
    else
        terminal_log=global_logger(ConsoleLogger(stderr, LogLevel(Info)))
    end

    log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"

    logger = TeeLogger(
        terminal_log,
        # Accept any messages with level >= Info
        MinLevelLogger(
            FileLogger("$output_path/logs/logfile$log_suffix.log"),
            Logging.Info
        ),
        # Accept any messages with level >= Debug
        MinLevelLogger(
            FileLogger("$output_path/logs/debug$log_suffix.log"),
            Logging.Debug,
        )
    )
    global_logger(logger)

end

# read JSON configuration file
function read_config(file_path::String)
"""
Read the JASON configuration file and parse it into a Dict
"""
    json_string = read(file_path, String)
    config = JSON.parse(json_string)
    return config
end

# process parsed arguments for the main function
function get_argparse()
"""
Parse the script arguments
"""
    settings = ArgParseSettings(prog="LEGEND ovbb Bayesian unbinned fit",
                            description="",
                            commands_are_required = true)
    @add_arg_table settings begin
        "--config", "-c"
            help = "path to config file"
            arg_type = String
            required = true
        "--index_toy", "-i"
            help = "index of sensitivity toy"
            arg_type = Int
            required = false
    end
    
    parse_args(settings)
    return parse_args(settings)
end

function main()
    
    # read parsed arguments
    @info "running using ",Base.Threads.nthreads()," threads"
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
   
    set_logger(config,output_path)

    # Call the analysis function from ZeroNuFit
    ZeroNuFit.run_analysis(config,output_path=output_path,toy_idx=nothing)

end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
