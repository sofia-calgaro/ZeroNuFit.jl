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


## Using already existing toys
Another way to run the code is present if, for instance, an user wants to use toy data generated according to one model but fit them with another model.
In this case, the path to the folder containing the already existing JSON files with toy data has to be provided together with the toy index:

```
julia sensitivity.jl -c config_fake_data.json -i N --path_to_toys path_to_your_toys
```

Below, an updated version of a bash file that can be used for retrieving multiple existing toy data and running sensitivity studies as multiple jobs on NERSC:

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

# set the directory path to toys
path_to_toys="output/fit_9_l200_1BI_new_data_coax_bkg_noS/sensitivity/fake_data" #"path/to/your/toys"
all_files=("$path_to_toys"/*.json)
full_paths=()
for file in "${all_files[@]}"; do
    if [[ -f "$file" ]]; then 
        full_paths+=("$file")
    fi
done
if [ ${#full_paths[@]} -eq 0 ]; then
    echo "The list of existing toy data is empty! Exit here."
    exit 1
else
    echo "You are going to run a fit over ${#full_paths[@]} number of already existing toys stored under $path_to_toys"
fi

# array to hold toy_idx
toy_indices=()

# Loop over available fake JSON toys
for path in "${full_paths[@]}"; do
    base_name="${path%.json}"
    number_str="${base_name##*fake_data}"  
    toy_idx=$((number_str))  

    toy_indices+=("$toy_idx") 
done
echo "List of toy indices: ${toy_indices[*]}"

# parallel execution - convert array to a space-separated list
module load parallel
module load julia
srun="srun -N 1"
parallel="parallel --delay 1 -j 128"
$srun $parallel "julia sensitivity.jl -c config/toy_9_l200_1BI_new_data_same_bkg_noS.json -i {1}" ::: "${toy_indices[@]}"

wait
```