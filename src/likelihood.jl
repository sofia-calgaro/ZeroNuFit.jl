using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets

using TypedTables
using Plots,LaTeXStrings
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
    b_name = part_k.bkg_name

    model_s_k = log(2) * N_A * part_k.exposure * (part_k.eff_tot + p.α * part_k.eff_tot_sigma) * (p.S*sig_units) / m_76
    model_b_k = deltaE * part_k.exposure * p[b_name]
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
    b_name = part_k.bkg_name
    model_s_k = log(2) * N_A * part_k.exposure * (part_k.eff_tot + p.α * part_k.eff_tot_sigma) * (p.S*sig_units) / m_76
    model_b_k = deltaE * part_k.exposure * p[b_name]
    model_tot_k = model_b_k + model_s_k
    
    # constrain λ not to be negative
    if model_tot_k < 0 
        λ = eps(0.)
    else
        λ = model_tot_k+eps(model_tot_k)
    end

    ll_value += logpdf(Poisson(λ), length(events_k))

  
    for evt_energy in events_k
        term1 = model_b_k / deltaE # background

        if (stat_only==true)
            term2 = model_s_k * pdf(Normal(Qbb + part_k.bias, part_k.fwhm/2.355), evt_energy) # signal (fixed nuisance)
        else
            term2 = model_s_k * pdf(Normal(Qbb + p.𝛥[idx_part_with_events], p.σ[idx_part_with_events]), evt_energy) # signal (free nuisance)
        end
        ll_value += log( (term1 + term2)+eps(term1+term2)) - log(model_tot_k+eps(model_tot_k)) 
    end

   
    
    return ll_value
end

# Tuple{Real, Real, Vector{Real}, Vector{Real}}
function build_likelihood_looping_partitions(partitions, events,part_event_index;stat_only=false,sqrt_prior=false,s_max=nothing)
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
    @debug part_event_index
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
            
            ## trick to include the prior 
            if (sqrt_prior)
                total_ll+=-log(2)-0.5*log(s_max)-0.5*log(p.S+eps(p.S))
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
    
    if config["signal"]["prior"] == "uniform" ||
        config["signal"]["prior"]=="sqrt"
        distrS = 0..uppS
    elseif config["signal"]["prior"]=="loguniform"
        distrS=LogUniform(0.01,uppS)
    
    else
        @error "distibution", config["signal"]["prior"], " is not yet defined"
        exit(-1)
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


    list_names = partitions.bkg_name
    unique_list=unique(list_names)

    bkg_par_names=[Symbol(name) for name in unique_list]

     
    distrS, distrB = get_signal_bkg_priors(config)
    distrB_multi=Dict(Symbol(bkg_par_name)=>distrB for bkg_par_name in bkg_par_names)

    
    
    if (stat_only==false)
        
        pretty_names =Dict(:S=>string("S [")*L"10^{-27}"*string("yr")*L"^{-1}"*string("]"),
        :α=>L"\alpha",
        :σ=>[],
        :𝛥=>[])

        res=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))
        bias=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))

        for (idx,part) in enumerate(partitions)
            
            if (part_event_index[idx]!=0)
                i_new = part_event_index[idx]
                res[i_new]=Truncated(Normal(part.fwhm/2.355,part.fwhm_sigma/2.355),0,Inf)
                bias[i_new] =Truncated(Normal(part.bias,part.bias_sigma),-Inf,Inf)
                long_name = string(part.experiment)*" "*string(part.part_name)*" "*part.detector
                append!(pretty_names[:σ],["Energy Resolution "*L"(\sigma)"*" "*long_name*" [keV]"])
                append!(pretty_names[:𝛥],["Energy Scale Bias "*L"(\Delta)"*" - "*long_name*" [keV]"])
    
            end
        end
        # get the minimum for α not to have negative values later on
        all_eff_tot = partitions.eff_tot
        all_eff_tot_sigma = partitions.eff_tot_sigma
        ratio = - all_eff_tot ./ all_eff_tot_sigma 
        α_min = maximum(ratio)
       

        # make some nice names for plotting
       
       

        for key in keys(distrB_multi)
            pretty_names[key]=string(key)*" [cts/keV/kg/yr]"
        end
        
        


        return distprod(S=distrS,;distrB_multi..., α=Truncated(Normal(0,1),α_min,Inf), σ=res, 𝛥=bias),pretty_names
        
    
    else 

        # simpler for a stat only fit
        pretty_names =Dict(:S=>string("S [")*L"10^{-27}"*string("yr")*L"^{-1}"*string("]"))
        for key in keys(distrB_multi)
            append!(pretty_names[key],[string(key)*" [cts/keV/kg/yr]"])
        end

        distprod(S=distS;distrB_multi...),
        Dict(:S=>L"S [10^{-27} \text{yr^{-1}}]",:B=>"B [cts/keV/kg/yr]")

    end
    
end


##############################################
##############################################
##############################################
function build_hd_prior(partitions,part_event_index;config,stat_only=false)
    """
    [experimental ] builds the priors for use in the fit with a Hierachical structure
    ----------
    Parameters
        - partitions:Table of the partition info
        - config: the Dict of the fit config
        - stat_only; a bool for whether systematic uncertatinties are considered on energy scale
    """
    
    
        list_names = partitions.bkg_name
        unique_list=unique(list_names)
    
        bkg_par_names=[Symbol(name) for name in unique_list]
    
         
        distrS, distrB = get_signal_bkg_priors(config)
        distrB_multi=Dict(Symbol(bkg_par_name)=>distrB for bkg_par_name in bkg_par_names)
    
        
        
        if (stat_only==false)
            
            pretty_names =Dict(:S=>string("S [")*L"10^{-27}"*string("yr")*L"^{-1}"*string("]"),
            :α=>L"\alpha",
            :σ=>[],
            :𝛥=>[])
    
            res=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))
            bias=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))
    
            for (idx,part) in enumerate(partitions)
                
                if (part_event_index[idx]!=0)
                    i_new = part_event_index[idx]
                    res[i_new]=Truncated(Normal(part.fwhm/2.355,part.fwhm_sigma/2.355),0,Inf)
                    bias[i_new] =Truncated(Normal(part.bias,part.bias_sigma),-Inf,Inf)
                    long_name = string(part.experiment)*" "*string(part.part_name)*" "*part.detector
                    append!(pretty_names[:σ],["Energy Resolution "*L"(\sigma)"*" "*long_name*" [keV]"])
                    append!(pretty_names[:𝛥],["Energy Scale Bias "*L"(\Delta)"*" - "*long_name*" [keV]"])
        
                end
            end
            # get the minimum for α not to have negative values later on
            all_eff_tot = partitions.eff_tot
            all_eff_tot_sigma = partitions.eff_tot_sigma
            ratio = - all_eff_tot ./ all_eff_tot_sigma 
            α_min = maximum(ratio)
           
    
            # make some nice names for plotting
           
           
    
            for key in keys(distrB_multi)
                pretty_names[key]=string(key)*" [cts/keV/kg/yr]"
            end
            
            dis_B = distprod
            hd = BAT.HierarchicalDistribution(
                    v -> begin 
                    dict = (; (key =>LogNormal(v.B,v.σB) for key in keys(distrB_multi))...)
                    BAT.NamedTupleDist(;dict...)
                    end,
                    BAT.NamedTupleDist(S=distrS,B=distrB,σB=distrB
                    , α=Truncated(Normal(0,1),α_min,Inf), σ=res, 𝛥=bias
                    )
            )          
            pretty_names[:B]="B [cts/keV/kg/yr]"
            pretty_names[:σB]=L"\sigma_B"*string("[cts/keV/kg/yr]")

    
            return hd,pretty_names
            
        
        else 
    
            # simpler for a stat only fit
            pretty_names =Dict(:S=>string("S [")*L"10^{-27}"*string("yr")*L"^{-1}"*string("]"))
            for key in keys(distrB_multi)
                append!(pretty_names[key],[string(key)*" [cts/keV/kg/yr]"])
            end
            pretty_names[:B]="B [cts/keV/kg/yr]"
            pretty_names[σB]=L"\sigma_B"*string("[cts/keV/kg/yr]")
            distprod(S=distS;distrB_multi...),
            Dict(:S=>L"S [10^{-27} \text{yr^{-1}}]",:B=>"B [cts/keV/kg/yr]")
    
        end
        
    end
    
    