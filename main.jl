### main.jl
### -> gets a config.json file in input for running the Bayesian unbinned fit
###

using Pkg
Pkg.activate(".") # Activate the environment

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
    
    ### TO DO: parse the config path to main()
    ### TO DO: if we save results somewhere (we can start from output/), it would be nice to give a name to saved files and/or stored them at specific paths - we can add this to the config file
    
    config_path = "config/config.json" 
    println("Reading configuration from: ", config_path)
    config = read_config(config_path)

    # Call the analysis function from ZeroNuFit
    ZeroNuFit.run_analysis(config)
end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
