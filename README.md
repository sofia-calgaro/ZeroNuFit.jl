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

```
{
    "debug":false,
    "partitions":["config/partitions_gerda_new.json","config/partitions_l200.json","config/partitions_mjd_new.json"],
    "events":    ["config/events_gerda.json","config/events_l200.json","config/events_mjd_new_part.json"],
    "output_path": "output/fit_mjd_l200_gerda_v2/",
    "overwrite":true,
    "bat_fit": {"nsteps": 1e6, "nchains": 6},
    "plot": {
            "fit_and_data": false,
            "bandfit_and_data": false,
            "scheme":"red",
            "alpha":0.3
        },
    "bkg_only": false,
    "signal": {"upper_bound":1000, "prior": "uniform"},
    "bkg": {"upper_bound":0.1, "prior": "uniform", "correlated": true},
    ...
}
```

where
- `"debug": true` if you want to display debug output on terminal;
- `"partitions"`: list of partitions JSON inputs; it takes one entry per experiment;
- `"events"`: list of events JSON inputs; it takes one entry per experiment;
- `"output_path"`: path where to store outputs (logs, plots, mcmc results);
- `"overwrite": true` if you want to overwrite a previous fit with same `output_path`; if set to `false` but no fits were previously performed (ie there are no outputs to overwrite), the code will save the output of this fit;
- `"bat_fit"`: settings for the BAT fit;
- `"plot"`: settings for plotting; `"fit_and_data": true` plots fit line over data (and CI bands if `"bandfit_and_data": true`); `"scheme":"red"` and `"alpha":0.3` are used for customizing output appearances;
- `"bkg_only": true` if we fit assuming no signal (S=0), false otherwise;
- `"signal"`: select `"upper_bound"` for the prior and the `"prior"` shape (`uniform`, `sqrt`, ...);
- `"bkg"`: select `"upper_bound"` for the prior and the `"prior"` shape (`uniform`, ...) and if you want to use a hierarchical model for correlations (`"correlated": true`).

Moreover, the config requires the following block for nuisance parameters, ie energy scale (=energy bias and resolution) and efficiency:
```
    {
    ...
    "nuisance": { 
        "energy_scale" : {
            "correlated": true,
            "fixed":     false
            },
         "efficiency" : {
            "correlated": true,
            "fixed": false
            }
    }
```

In particular, you can set `"correlated": true` if you want to use one variable to correlate the nuisance parameters (eg to speed up the computation times), and `"fixed": false` if you want to include a prior for nuisance parameters (otherwise these parameters they will be fixed to their partition value and not constrained).
 
If a variable is correlated (either `energy_scale` or `efficency`), the code will search for a field in the `fit_groups` block of the partitions JSON file to use a correlated variable per each fit group. 
In particular, the field has to be specified as:
- `"efficiency_group_name": "..."`
- `"energy_scale_group_name": "..."`

 > [!NOTE] 
 > If the key doesn't exist, this defaults to "all"
 
Parameters are then added to the model called `αr_$name` (for resolution), `αe_$name` for efficiency and `αb_$name` for bias.
 
 > [!WARNING]
 > The $\alpha$ parameter names default to `_all`, if you want one different per experiment this must be explicitly specified in the fit groups entry


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
To convert to this format from the standard GERDA and LEGEND files, there are tools available in https://github.com/tdixon97/legend-0vbb-config.

It is possible to supply a list of partition and event files in this case the list of fit groups and events are concatenated.

> [!WARNING]  
> If multiple files are provided `fit_group` must still be unique.


## Sensitivity studies
Another module is present for running sensitivity studies. This can be run as

```
julia sensitivity.jl -c config/config_fake_data.json -i N
```

where `N` is an integer number corresponding to the toy index.
The command can be run in an external bash script for looping over this index.

The input config file has the following entries:

```
{
    "path_to_fit": "output/fit_alpha_high_stat_true_TOBY4_v4/",
    "best_fit": false,
    "seed": null
}

```

where
- `"path_to_fit"` is the path to the already performed fit over real data;
- `"best_fit": true` if we want to fix the paramaters to the best fit;
- `"seed": null` if we want a random seed when generating fake data, otherwise you can fix it to an Int value.

Any information about the signal being included or not in the fit of real data, was saved and retrieved from the output JSON file with results.