### main.jl
### -> gets a config.json file in input for running the Bayesian unbinned fit

using JSON

# load the script to run the analysis
include("src/ZeroNuFit.jl")
using .ZeroNuFit


# Function to read JSON configuration file
function read_config(file_path::String)
    json_string = read(file_path, String)
    config = JSON.parse(json_string)
    return config
end

function main()

    # Check if the config argument was provided
    #if !haskey(args, :config)
    #    println("Error: Configuration file not provided. Use -c option.")
    #    return
    #end

    config_path = "config/config.json" # how to argparse from the outside?
    println("Reading configuration from: ", config_path)

    config = read_config(config_path)

    # Call the analysis function from ZeroNuFit
    ZeroNuFit.run_analysis(config)
end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
