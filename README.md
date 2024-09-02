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
