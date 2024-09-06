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

## Partition and events files
The takes inputs in JSON format, two files are needed a "partitions file" giving information on the independent spectra to be used in the fit/likelihood, this is set by the "partitions" key in the config file. This provides all the information neccesary to define the fit model.

The file consists of a file of independent spectra to include in the fit (for example channels or partitions). A partition is defined uniquely by a range of time-stamps, a detector name and an experiment name. 

> [!NOTE]
> In principle the 'detector' does not need to be a single detector but can be a label for any groups of detectors. This allows to make fits where all detectors are fit together.

The partitions are grouped into `fit_groups` these are sets of partitions which are treated with the same background fit model and range.
In the partitions file the user must provide the information on the fit groups and partitions, (organised by fit group). This file must be provided as a JSON file, this allows a full customisation of the fit.

This JSON file has a nested structure with two subdictonaries, the first with key "fit_groups", describing the groupings in the fit, and the second "partitions" giving a list of partitions for each fit group.
An example is shown below.
```
{
"fit_groups":{
                "group_one":{
                            "range":[[1930,2099],[2109,2114],[2124,2190]],
                            "model":"uniform",
                            "bkg_name":"low_bkg"
                            }

},
"partitions": {
                "group_one":[

                            {  
                                "experiment": "LEGEND",
                                "detector": "DET_0",
                                "start_ts": 1704950367,
                                "end_ts": 1708271505,
                                "eff_tot": 0.6,
                                "eff_tot_sigma": 0.1,
                                "fwhm": 3,
                                "fwhm_sigma": 1,
                                "exposure": 1,
                                "bias": 0.2,
                                "bias_sigma": 0.1
                            }, ...
                            ],
                "group_two":...
            },


}
            
```
in future we will also add the possibility to customize further the fit. Currently it implements a fit to the energy spectrum with a uniform background.
They key `bkg_name` is used to set the name of the background parameter for this group, note that several groups can be fitted with the same background parameter, this enables quick modifcation of the fit.

In addition, it is neccesary to provide an 'event' file describing the events observed in the data, the path to this file is specified by the 'events' key in the config. Again this is a JSON file consisting of a list of observed events of the form.
 
```
    {       "experiment":"LEGEND",
            "energy": 2069.420,
            "timestamp": 1755109448,
            "detector": "DET_0"
        },
```
The timestamp and detector are used to extract which partition this event corresponds to.
To convert to this format from the standard GERDA and LEGEND files , there is a notebook called `make_configs.ipynb' containing the neccesary functions.

It is possible to supply a list of partition and event files in this case the list of fit groups and events are concatenated.

> [!WARNING]  
> If multiple files are provided `fit_group` must still be unique.