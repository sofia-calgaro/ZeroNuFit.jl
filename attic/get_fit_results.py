"""
Script to load main output results from a fit
Main Authors: Sofia Calgaro, Toby Dixon
"""
import os,json
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Script to load fit results"
    )
    parser.add_argument(
        "--fit", "--f",
        type=str,
        help="Name of the fit output folder set initially in the config file (where results are stored)",
    )
    parser.add_argument(
        "--path", "--p",
        type=str,
        default="output",
        help="Output path: 'output' is default",
    )
    parser.add_argument(
        "--nuisance", "--n",
        default="False",
        type=str,
        help="Plot nuisance parameters (i.e. efficiencies, resolutions, biases) output results? 'False' is deafult (case sensitive)",
    )
    parser.add_argument(
        "--global_mode", "--gm",
        default="refined_global_modes",
        type=str,
        help="Global mode type: 'refined_global_modes' (default) or 'global_modes'",
    )
    args = parser.parse_args()
    nuisance = False if args.nuisance == "False" else True
    gm_type = args.global_mode
    fit = args.fit
    path = args.path
    
    print(f"You are inspecting '{fit}'")

    # path were you can find the fit results
    full_path = os.path.join(path, fit, "mcmc_files/fit_results.json")
    json_file = json.load(open(full_path))

    quantiles = json_file['quantile90']
    gms = json_file[gm_type]
    c68 = json_file['ci_68']

    print("*****************")
    print("SIGNAL UPPER LIMIT")
    if 'S' in quantiles.keys():
        print(f" T1/2 > {1/(quantiles['S']*1e-27)} x 10^27 yr")
        print(f" S < {quantiles['S']} x 10^-27 yr^-1")

    print("")
    print("*****************")
    print("SIGNAL GLOBAL MODE")
    if 'S' in quantiles.keys():
        if gms['S'] != 0:
            t12 = 1/(gms['S']*1e-27)
            print(f" T1/2 = {t12} x 10^27 yr")
            print(f" S = {gms['S']} x 10^-27 yr^-1")

    print("")
    print("*****************")
    print("BACKGROUND PARAMETERS")
    for par in gms.keys():
        if "B" not in par: continue
        marg_par = json_file['marginalized_modes'][par]*1e4
        gm_par = json_file['refined_global_modes'][par]*1e4
        c68low = c68[par][0]['left']*1e4
        c68high = c68[par][0]['right']*1e4
        print("")
        print(par)
        print("Marginalized mode:", marg_par, "1e-4 ckky")
        print("Global mode:", gm_par, "1e-4 ckky")
        print("Smallest 68% CI:", f"[{c68low}, {c68high}]", "1e-4 ckky")

    for par in gms.keys():
        if "ph" not in par: continue
        marg_par = json_file['marginalized_modes'][par]*1e4
        gm_par = json_file['refined_global_modes'][par]*1e4
        c68low = c68[par][0]['left']*1e4
        c68high = c68[par][0]['right']*1e4
        print("")
        print(par)
        print("Marginalized mode:", marg_par, "1e-4 ckky")
        print("Global mode:", gm_par, "1e-4 ckky")
        print("Smallest 68% CI:", f"[{c68low}, {c68high}]", "1e-4 ckky")

    if nuisance:
        print("")
        print("*****************")
        print("NUISANCE PARAMETERS")
        for par in gms.keys():
            if "S" in par or "ph" in par or "B" in par: continue

            # the parameter has 1 entry
            if not isinstance(gms[par], list):
                marg_par = json_file['marginalized_modes'][par]
                gm_par = json_file['refined_global_modes'][par]
                c68low = c68[par][0]['left']
                c68high = c68[par][0]['right']
                print("")
                print(par)
                print("Marginalized mode:", marg_par)
                print("Global mode:", gm_par)
                print("Smallest 68% CI:", f"[{c68low}, {c68high}]")
            # the parameter is a list of entries
            else:
                for idx,entry in enumerate(gms[par]):
                    marg_par = json_file['marginalized_modes'][par]
                    gm_par = json_file['refined_global_modes'][par]
                    c68low = c68[par][idx][0]['left']
                    c68high = c68[par][idx][0]['right']
                    print("")
                    print(par, "-", idx)
                    print("Marginalized mode:", marg_par)
                    print("Global mode:", gm_par)
                    print("Smallest 68% CI:", f"[{c68low}, {c68high}]")

if __name__ == "__main__":
    main()