using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets
using TypedTables
using Plots
using Cuba

# define some constants -TODO add to config
Qbb = 2039.06 # keV
N_A = 6.022E23
m_76 = 75.6E-3 # kg/mol
deltaE = 240 # keV
sig_units =1e-27 # signal is in units of this


function build_likelihood_zero_obs_evts(part_k, p;stat_only=false)
"""
Function to calculate the partial likelihood for a partition with 0 events
    
"""

    ll_value = 0
    
    model_s_k = log(2) * N_A * part_k.exposure * (part_k.eff_tot + p.α * part_k.eff_tot_sigma) * (p.S*sig_units) / m_76
    model_b_k = deltaE * part_k.exposure * p.B
    model_tot_k = model_b_k + model_s_k

    ll_value += -(model_tot_k+eps(model_tot_k)) 
    
    return ll_value
end

function build_likelihood_per_partition(idx_k, idx_part_with_events,part_k, events_k, p;stat_only=false)
"""
Function which computes the partial likelihood for a single data partiton
free parameters: signal (S), background (B), energy bias (biask) and resolution per partition (resk)
"""

    ll_value = 0
    
    model_s_k = log(2) * N_A * part_k.exposure * (part_k.eff_tot + p.α * part_k.eff_tot_sigma) * (p.S*sig_units) / m_76
    model_b_k = deltaE * part_k.exposure * p.B
    model_tot_k = model_b_k + model_s_k
    
    # constrain λ not to be negative
    if model_tot_k < 0 
        λ = eps(0.)
    else
        λ = model_tot_k+eps(model_tot_k)
    end

    ll_value += logpdf(Poisson(λ), length(events_k))
    
    for i in events_k
        for evt_energy in events_k

            term1 = model_b_k / deltaE # background

            if (stat_only==true)
                term2 = model_s_k * pdf(Normal(Qbb + part_k.bias, part_k.fwhm/2.355), evt_energy) # signal (fixed nuisance)
            else
                term2 = model_s_k * pdf(Normal(Qbb + p.𝛥[idx_part_with_events], p.σ[idx_part_with_events]), evt_energy) # signal (free nuisance)
            end
            ll_value += log( (term1 + term2)+eps(term1+term2)) - log(model_tot_k+eps(model_tot_k)) 
        end

    end
    
    return ll_value
end

# Tuple{Real, Real, Vector{Real}, Vector{Real}}
function build_likelihood_looping_partitions(partitions, events,part_event_index;stat_only=false)
"""
Function which creates the likelihood function for the fit (looping over partitions)
Parameters:
-----------
    -partitions: Table - partitions input file
    -events: Array      - list of events in each partitions (with their energy)
    -stat_only:bool     -whether the fit includes only parameters of interest
Returns:
--------
    DensityInterface.logfuncdensity - the likelihood function
"""

    return DensityInterface.logfuncdensity( function(p)
            total_ll = 0.0
            
            for (idx_k, part_k) in enumerate(partitions)
                
                if part_event_index[idx_k]!=0
                    idx_k_with_events=part_event_index[idx_k]
                    total_ll += build_likelihood_per_partition(idx_k,part_event_index[idx_k], part_k, events[idx_k], p, stat_only=stat_only)
                else
                    # no events are there for a given partition
                    total_ll += build_likelihood_zero_obs_evts(part_k, p)
                end
            end
            return total_ll
        end
    )
end


##############################################
##############################################
##############################################
function get_signal_bkg_priors(config)
"""
Defines specific priors for signal and background contributions
Parameters
    - config: the Dict of the fit config
"""
    
    uppS = config["signal"]["upper_bound"]
    uppB = config["bkg"]["upper_bound"]
    
    if config["signal"]["prior"] == "uniform"
        distrS = 0..uppS
    end
    if config["bkg"]["prior"] == "uniform"
        distrB = 0..uppB
    end
    
    return distrS, distrB
end

##############################################
##############################################
##############################################
function build_prior(partitions,part_event_index;config,stat_only=false)
"""
Builds the priors for use in the fit
----------
Parameters
    - partitions:Table of the partition info
    - config: the Dict of the fit config
    - stat_only; a bool for whether systematic uncertatinties are considered on energy scale
"""
    res=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))
    bias=Vector{Normal{Float64}}(undef,maximum(part_event_index))

    for (idx,part) in enumerate(partitions)
        
        if (part_event_index[idx]!=0)
            i_new = part_event_index[idx]
            res[i_new]=Truncated(Normal(part.fwhm/2.355,part.fwhm_sigma/2.355),0,Inf)
            bias[i_new] =Normal(part.bias,part.bias_sigma)
        end
    end
    if (stat_only==false)
        # get the minimum for α not to have negative values later on
        all_eff_tot = partitions.eff_tot
        all_eff_tot_sigma = partitions.eff_tot_sigma
        ratio = - all_eff_tot ./ all_eff_tot_sigma 
        α_min = maximum(ratio)
        
        distrS, distrB = get_signal_bkg_priors(config)
        return distprod(S=distrS,B=distrB, α=Truncated(Normal(0,1),α_min,Inf), σ=res, 𝛥=bias)
    
    else 
        distprod(S=0..config["signal"]["upper_bound"],B=0..config["bkg"]["upper_bound"])
    end
    
end

