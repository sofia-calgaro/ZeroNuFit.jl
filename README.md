A Bayesian unbinned fit based on [BAT.jl](https://github.com/bat/BAT.jl)

## How to run the code
Run the following command by specifying the path to the configuration file used for settings:

```
julia main.jl -c config/config.json
```

## Julia Project enviroments
To run the code in a virtual enviroment you can use the following.
```
julia
] 
activate .
instantiate

```
Now you can run the script inside this enviroment with:

```
julia main.jl --project=. -c config/config.json
```
alternatively the code can be run within the Legend container, see https://github.com/legend-exp/legend-julia-tutorial for more details.


## Config file
Before running the code, set the input config.json file with following entries:

```bash
{
    "debug":false, // true if you want to display debug output on terminal
    "partitions":["config/partitions_gerda.json"], // include partitions inputs -> one entry per experiment
    "events":    ["config/events_gerda.json"], // include events inputs -> one entry per experiment
    "meta_path": "/global/cfs/cdirs/m2676/data/lngs/l200/public/prodenv/prod-blind/ref-v2.1.2/inputs", // path to metadata
    "output_path": "output/test_fit/", // path for storing outputs (logs, plots, mcmc results)
    "bat_fit": {"nsteps": 1e3, "nchains": 4}, // some settings for running the BAT fit
    "upper_signal":1e-25, // upper bound for signal (free parameter)
    "upper_bkg":1e-3, // upper bound for background (free parameter)
    "stat_only":true // false if we use pull terms for additional nuisance parameters (res, bias, eff.)
}
```
