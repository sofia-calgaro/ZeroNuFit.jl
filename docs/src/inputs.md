# Partitions and events

The fit takes in inputs two files in JSON format (for a full customization of the fit), which paths have to be specified in the `config.json` file.

Table of contents:

```@contents
Pages = ["inputs.md"]
Depth = 3
```


## Partitions file
The partitions file gives information on the independent spectra to be used in the fit/likelihood, this is set by the "partitions" key in the config file. 
This provides all the information neccessary to define the fit model.

The file consists of a file of independent spectra to include in the fit (for example channels or partitions). 
A partition is defined uniquely by a range of time-stamps, a detector name and an experiment name. 

!!! note

    In principle the 'detector' does not need to be a single detector (e.g. `C000RG1`) but can be a label for any groups of detectors (e.g. `COAX_dets`). 
    This allows to make fits where all detectors are fit together.

The partitions are grouped into `fit_groups`: these are sets of partitions which are treated with the same background/signal fit model and range.
In the partitions file, the user must provide the information on the fit groups and partitions (organized by fit group). 

This JSON file has a nested structure with two subdictionaries, the first with key `"fit_groups"`, describing the groupings in the fit, and the second `"partitions"` giving a list of partitions for each fit group.
An example is shown below.
```
{
"fit_groups":{
                "group_one":{
                            "range":[[1930,2099],[2109,2114],[2124,2190]],
                            "model":"uniform",
                            "bkg_name":"low_bkg"
                            "signal_name":"gaussian_plus_lowEtail"
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
                                "width": 3,
                                "width_sigma": 1,
                                "exposure": 1,
                                "bias": 0.2,
                                "bias_sigma": 0.1
                            }, ...
                            ],
                "group_two":...
            },


}
            
```
They key `bkg_name` is used to set the name of the background parameter for this group.
Note that several groups can be fitted with the same background parameter, this enables quick modification of the fit.

The key `"model":"uniform"` is used to set the background model to uniform as default.
For different background model shapes, additional information are necessary and these can be specified in the `config.json` (see the "Configuration file" documentation).

!!! warning

    The background shape is set to global for all partitions, differently from the signal shape (see below) that can be specified differently for each fit group.

Notice it is possible to specify the signal shape for each fit group. 
The available options at the moment are `"signal_name":"gaussian_plus_lowEtail"` or `"signal_name":"gaussian"` (default).
If the key is omitted, the default Gaussian signal shape will be adopted.
Notice that if you want to use the `"signal_name":"gaussian_plus_lowEtail"` option (eg for MAJORANA DEMONSTRATOR data), you need to provide additional input signal shape parameters (`frac`, `sigma`, `tau`). 

!!! note

    MAJORANA DEMONSTRATOR (MJD) input data are taken from ["I. J. Arnquist et al., Final Result of the Majorana Demonstrator’s Search for Neutrinoless Double- β Decay in Ge 76, PRL 130, 062501 (2023)"](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.062501). 
    Here, you can find the input data under the supplemental materials: `supp_analysis_parameters.txt` and `supp_event_list.txt`.



!!! warning

    Notice that the `width` parameter assumes different meanings depending on the chosen signal shape:
    - if `"signal_name": "gaussian_plus_lowEtail"`, `width` is the fractional uncertainty on the FWHM of the peak at 2039 keV (it rescales the energy resolution of the peak)
    - if `"signal_name": "gaussian"`, `width` is the energy resolution expressed in standard deviations (NOT in FWHM)


## Events file
In addition, it is neccessary to provide an 'event' file describing the events observed in the data, the path to this file is specified by the 'events' key in the config. Again this is a JSON file consisting of a list of observed events of the form.
 
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

!!! warning

    If multiple files are provided `fit_group` must still be unique.

!!! note

    Two example scripts of how to build the events and partitions files are present under `attic`, see `make_input_events.jl` and `make_input_partitions.jl`, respectively. By adjusting the functions to your needs, you will be able to easily populate the content of JSON files you need to load in the config file before running your fit. 
    