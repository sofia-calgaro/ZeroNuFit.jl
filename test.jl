
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
using DensityInterface
using BenchmarkTools
using Debugger


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
    part_event_index,events,partitions= get_partitions_events(config)
    prior,likelihood,posterior,par_names=get_stat_blocks(partitions,events,part_event_index,config=config,bkg_only=false)

    /legend-0vbb-config/partitions_gerda_pI.json","../legend-0vbb-config/partitions_gerda_fix_bias.json", "../legend-0vbb-config/partitions_l200.json", "../legend-0vbb-config/partitions_mjd_new.json
    pars = (

    # global pars
    reso = fill(-0.01, 1)
    append!(reso,fill(1.8, 1))
    append!(reso,)
    S=0.5,
    αe_all=0.01,

    B=1E-2,
    σ=[fill(-0.01, 46),]
    gerdaI_golden_Δk=fill(-0.01, 46),
    gerdaI_golden_σk=fill(1.8, 46),

    gerdaI_silver_B=1E-3,
    gerdaI_silver_Δk=fill(0.01, 10),
    gerdaI_silver_σk=fill(1.8, 10),

    gerdaI_bege_B=5E-4,
    gerdaI_bege_Δk=fill(0.02, 3),
    gerdaI_bege_σk=fill(1.1, 3),

    gerdaI_extra_B=3E-4,
    gerdaI_extra_Δk=fill(-0.02, 2),
    gerdaI_extra_σk=fill(1.7, 2),

    gerdaII_B=4E-4,
    gerdaII_Δk=fill(0.01, 13),
    gerdaII_σk=fill(1.2, 13),

    majorana_DS0_B=2E-2,
    majorana_DS0_Δk=fill(0.03, 7),
    majorana_DS0_σk=fill(2.5, 7),

    majorana_mod1_B=1E-2,
    majorana_mod1_Δk=fill(-0.03, 75),
    majorana_mod1_σk=fill(1.1, 75),

    majorana_mod2_B=4E-2,
    majorana_mod2_Δk=fill(0.05, 14),
    majorana_mod2_σk=fill(1.1, 14),

    legend200_B=2E-4,
    legend200_Δk=fill(0.05, 7),
    legend200_σk=fill(1.1, 7),
)

end







if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
