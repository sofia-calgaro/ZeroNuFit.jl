# Generating toys
Another module is present for generating toys and running sensitivity studies. This can be run as

```
julia sensitivity.jl -c config_fake_data.json -i N
```

where `N` is an integer number corresponding to the toy index.
The command can be run in an external bash script for looping over this index.

The input config file (`config_fake_data.json`) has the following entries:

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

Below, we show an example of bash file used for running sensitivity studies as multiple jobs on NERSC:

```bash
#!/bin/bash                                                                                                                                                 
#SBATCH -q regular                                                                                                                                       
#SBATCH --constraint=cpu                                                                                                                                    
#SBATCH -t 48:00:00
#SBATCH -J sens_test                                                                                                                                         

#SBATCH --mail-user=<your_email>
#SBATCH --mail-type=ALL                                                                                                                                     
#SBATCH --output output_path/parallel.log                                                     
#SBATCH --error output_path/parallel.err  

#SBATCH  --image=legendexp/legend-base:latest               

module load parallel
module load julia
srun="srun -N 1"
parallel="parallel --delay 1 -j 128"

# run parallel jobs
$srun  $parallel "julia sensitivity.jl -c config_fake_data.json -i {1}" ::: {1..10000} &

wait
```