using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets

using TypedTables
using Plots,LaTeXStrings
using Cuba
using OrderedCollections

function get_mu_s_b(p::NamedTuple,part_k::NamedTuple,idx_part_with_events::Int;nuis_correlated::Bool=true,nuis_prior::Bool=false,bkg_only::Bool=false)
    """
    Get the expected number of signal and background counts in a partition
    """
    N_A = 6.022E23
    m_76 = 75.92E-3 # kg/mol
    deltaE = 240 # keV
    sig_units =1e-27 # signal is in units of this
    b_name = part_k.bkg_name

    model_s_k = 0
    if bkg_only==false
        if nuis_correlated == true
            model_s_k = log(2) * N_A * part_k.exposure * (part_k.eff_tot + p.α * part_k.eff_tot_sigma) * (p.S*sig_units) / m_76
        # we remove alpha and uncertainties
        else
            if (idx_part_with_events!=0 || nuis_prior==true)
                model_s_k = log(2) * N_A * part_k.exposure * p.ε[idx_part_with_events] * (p.S*sig_units) / m_76
            else
                model_s_k = log(2) * N_A * part_k.exposure * part_k.eff_tot * (p.S*sig_units) / m_76
            end
        end
    end
    model_b_k = deltaE * part_k.exposure * p[b_name]

    return model_s_k,model_b_k
end

function build_likelihood_zero_obs_evts(part_k::NamedTuple, p::NamedTuple;nuis_prior::Bool=false, nuis_correlated::Bool=true, bkg_only::Bool=false)
"""
Function to calculate the partial likelihood for a partition with 0 events
    
"""

    ll_value = 0
    model_s_k,model_b_k = get_mu_s_b(p,part_k,0,nuis_correlated=nuis_correlated,nuis_prior=nuis_prior,bkg_only=bkg_only)
    model_tot_k = model_b_k + model_s_k

    ll_value += -(model_tot_k+eps(model_tot_k)) 
    
    return ll_value
end

function build_likelihood_per_partition(idx_k::Int, idx_part_with_events::Int,part_k::NamedTuple, events_k::Vector{Float64}, p::NamedTuple;
                                        nuis_prior::Bool=false, nuis_correlated::Bool=true, bkg_only::Bool=false)
"""
Function which computes the partial likelihood for a single data partiton
free parameters: signal (S), background (B), energy bias (biask) and resolution per partition (resk)
"""
    Qbb = 2039.06 # keV
    deltaE = 240 # keV

    ll_value = 0

    model_s_k,model_b_k =   get_mu_s_b(p,part_k,idx_part_with_events,nuis_correlated=nuis_correlated,nuis_prior=nuis_prior,bkg_only=bkg_only)
   
    model_tot_k = model_b_k + model_s_k
    
    # constrain λ not to be negative
    if model_tot_k < 0 
        λ = eps(0.)
    else
        λ = model_tot_k+eps(model_tot_k)
    end

    ll_value += logpdf(Poisson(λ), length(events_k))

    if bkg_only==false
        for evt_energy in events_k
            term1 = model_b_k / deltaE # background

            if (nuis_prior==false)
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
function build_likelihood_looping_partitions(partitions::TypedTables.Table, events::Array{Vector{Float64}},part_event_index::Vector{Int};
                                            nuis_prior=true,nuis_correlated=true,sqrt_prior=false,s_max=nothing,bkg_only=false)
"""
Function which creates the likelihood function for the fit (looping over partitions)
Parameters:
-----------
    -partitions: Table - partitions input file
    -events: Array      - list of events in each partitions (with their energy)
    -nuis_prior:bool     - true if we want to include priors for nuisance parameters (bias, res, eff)
Returns:
--------
    DensityInterface.logfuncdensity - the likelihood function
"""
    @debug part_event_index
    return DensityInterface.logfuncdensity( function(p::NamedTuple)
            total_ll = 0.0
            
             for (idx_k, part_k) in enumerate(partitions)
                
                if part_event_index[idx_k]!=0
                    idx_k_with_events=part_event_index[idx_k]
                    total_ll += build_likelihood_per_partition(idx_k,part_event_index[idx_k], part_k, events[idx_k], p, nuis_prior=nuis_prior, nuis_correlated=nuis_correlated, bkg_only=bkg_only)
                else
                    # no events are there for a given partition
                    total_ll += build_likelihood_zero_obs_evts(part_k, p,nuis_prior=nuis_prior, nuis_correlated=nuis_correlated, bkg_only=bkg_only)
                end
            end
            
            ## trick to include thes sqrt prior 
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
function generate_data(samples::BAT.DensitySampleVector,partitions::TypedTables.Table,part_event_index::Vector{Int};
    best_fit::Bool=false,nuis_prior=true,bkg_only=false,seed=nothing,nuis_correlated=true)
"""
Generates data from a posterior distribution.
This is based on the posterior predictive distributions. 
Given a model with some parameters `theta_i`, the posterior predictive distribution,
or the distribution of data generated according to the posterior distribution of theta
and the likelihood is:
```math
p(y|D) =int p(y|theta)p(theta|D)dtheta
```
Or in terms of sampling we first draw samples of `theta` from the posterior and then generate,
datasets based on the likelihood.
We also give the options to fix the posterior distribution to the best fit,
which is equivalent to the standard sampling methods.

Parameters
----------
    - samples::DensitySamplesVector the samples of a past fit or a NamedTuple of best fit
    - partitions::Table of the partition info
    - part_event_index: index for the parameters for partitions with events
Keyword arguments
-----------------
    - best_fit::Bool where to fix the paramaters to the best fit
    - nuis_prior::Bool whether only statistical parameters were included in the posterior
    - bkg_only::Bool where the fit was without signal,
    - seed::Int random seed
Returns
    OrderedDict of the data
"""
    Qbb = 2039.06 # keV

    # seed the seed
    output=OrderedDict("events"=>[])
    if (seed==nothing)
        Random.seed!(Int(round(10000*(rand()))))
    else
        Random.seed!(seed)
    end
    
    
    if (samples isa NamedTuple)
        p= samples
    else
        distribution = Categorical(samples.weight/sum(samples.weight)) # Generate a random index based on the weights
        random_index = rand(distribution)
        p = samples.v[random_index]
    end
    # create the array to fill
    events=[]

    for (idx_k, part_k) in enumerate(partitions)

        b_name = part_k.bkg_name           
        idx_part_with_events=part_event_index[idx_k]
        model_s_k,model_b_k = get_mu_s_b(p,part_k,idx_part_with_events,nuis_correlated=nuis_correlated,nuis_prior=nuis_prior,bkg_only=bkg_only)

        n_s = rand(Poisson(model_s_k))
        n_b = rand(Poisson(model_b_k))
        events =generate_disjoint_uniform_samples(n_b)
        if (bkg_only == false)
            for i in 1:n_s
                if (nuis_prior==false || idx_part_with_events==0)
                    append!(events,rand(Normal(Qbb + part_k.bias, part_k.fwhm/2.355)))
                else    
                    append!(events,rand(Normal(Qbb + p.𝛥[idx_part_with_events], p.σ[idx_part_with_events])))
                end
        
            end
        end
        times = rand(Uniform(part_k.start_ts,part_k.end_ts),length(events))
        times=[Int(round(t)) for t in times]
        
        if (length(times)>0)
            for (t,e) in zip(times,events)
                append!(output["events"],[OrderedDict("timestamp"=>t,"experiment"=>part_k.experiment,
                "detector"=>part_k.detector,"energy"=>e
                )])
            end

        end

    end
    display(output["events"])
    
    return output

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
function build_prior(partitions,part_event_index;config,nuis_prior=true,bkg_only=false)
"""
Builds the priors for use in the fit
----------
Parameters
    - partitions:Table of the partition info
    - config: the Dict of the fit config
    - nuis_prior; true if we want to include priors for nuisance parameters (bias, res, eff)
"""


    list_names = partitions.bkg_name
    unique_list=unique(list_names)
    nuis_correlated = config["nuisances"]["correlated"]

    bkg_par_names=[Symbol(name) for name in unique_list]
     
    distrS, distrB = get_signal_bkg_priors(config)
    distrB_multi=Dict(Symbol(bkg_par_name)=>distrB for bkg_par_name in bkg_par_names)

    pretty_names =Dict(:S=>string("S [")*L"10^{-27}"*string("yr")*L"^{-1}"*string("]"),
        :α=>L"\alpha",
        :ε=>[],
        :σ=>[],
        :𝛥=>[])
    
    if (nuis_prior==true)
        
        # model efficiencies with an alpha parameter (if set to True)
        if nuis_correlated == true
            
            @info "...CORRELATED EFF IS TRUE!"

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

            # make some nice names for plotting
            for key in keys(distrB_multi)
                pretty_names[key]=string(key)*" [cts/keV/kg/yr]"
            end

            # get the minimum for α not to have negative values later on
            all_eff_tot = partitions.eff_tot
            all_eff_tot_sigma = partitions.eff_tot_sigma
            ratio = - all_eff_tot ./ all_eff_tot_sigma 
            α_min = maximum(ratio)
            
            if bkg_only==false
                return distprod(S=distrS,;distrB_multi..., α=Truncated(Normal(0,1),α_min,Inf), σ=res, 𝛥=bias),pretty_names
            else
                return distprod(;distrB_multi..., α=Truncated(Normal(0,1),α_min,Inf), σ=res, 𝛥=bias),pretty_names
            end
            
        else
            
            @info "...CORRELATED EFF IS FALSE!"

            res=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))
            bias=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))
            eff=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))

            for (idx,part) in enumerate(partitions)

                if (part_event_index[idx]!=0)
                    i_new = part_event_index[idx]
                    res[i_new]=Truncated(Normal(part.fwhm/2.355,part.fwhm_sigma/2.355),0,Inf)
                    bias[i_new] =Truncated(Normal(part.bias,part.bias_sigma),-Inf,Inf)
                    eff[i_new] =Truncated(Normal(part.eff_tot,part.eff_tot_sigma),0,1)
                    long_name = string(part.experiment)*" "*string(part.part_name)*" "*part.detector
                    append!(pretty_names[:σ],["Energy Resolution "*L"(\sigma)"*" "*long_name*" [keV]"])
                    append!(pretty_names[:𝛥],["Energy Scale Bias "*L"(\Delta)"*" - "*long_name*" [keV]"])
                    append!(pretty_names[:ε],["Efficiency "*L"(\varepsilon)"*" - "*long_name*""])
                end
            end

            # make some nice names for plotting
            for key in keys(distrB_multi)
                pretty_names[key]=string(key)*" [cts/keV/kg/yr]"
            end
            
            if bkg_only==false
                return distprod(S=distrS,;distrB_multi..., ε=eff, σ=res, 𝛥=bias),pretty_names
            else
                return distprod(;distrB_multi..., ε=eff, σ=res, 𝛥=bias),pretty_names
            end
        end
        
    
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
function build_hd_prior(partitions,part_event_index;config,nuis_prior=true,bkg_only=false)
"""
[experimental ] builds the priors for use in the fit with a Hierachical structure
----------
Parameters
    - partitions:Table of the partition info
    - config: the Dict of the fit config
    - nuis_prior; true if we want to include priors for nuisance parameters (bias, res, eff)
"""

    list_names = partitions.bkg_name
    unique_list=unique(list_names)
    nuis_correlated = config["nuisances"]["correlated"]

    bkg_par_names=[Symbol(name) for name in unique_list]

    distrS, distrB = get_signal_bkg_priors(config)
    distrB_multi=Dict(Symbol(bkg_par_name)=>distrB for bkg_par_name in bkg_par_names)

    pretty_names =Dict(:S=>string("S [")*L"10^{-27}"*string("yr")*L"^{-1}"*string("]"),
        :α=>L"\alpha",
        :ε=>[],
        :σ=>[],
        :𝛥=>[])

    if (nuis_prior==true)
        
        if nuis_correlated == true

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
            if bkg_only==false
                hd = BAT.HierarchicalDistribution(
                        v -> begin 
                        dict = (; (key =>LogNormal(log(v.B)-0.5*v.σB*v.σB,v.σB) for key in keys(distrB_multi))...)
                        BAT.NamedTupleDist(;dict...)
                        end,
                        BAT.NamedTupleDist(S=distrS,B=distrB,σB=0..1
                        , α=Truncated(Normal(0,1),α_min,Inf), σ=res, 𝛥=bias
                        )
                ) 
            else
                hd = BAT.HierarchicalDistribution(
                        v -> begin 
                        dict = (; (key =>LogNormal(log(v.B)-0.5*v.σB*v.σB,v.σB) for key in keys(distrB_multi))...)
                        BAT.NamedTupleDist(;dict...)
                        end,
                        BAT.NamedTupleDist(B=distrB,σB=0..1
                        , α=Truncated(Normal(0,1),α_min,Inf), σ=res, 𝛥=bias
                        )
                ) 
            end
            pretty_names[:B]="B [cts/keV/kg/yr]"
            pretty_names[:σB]=L"\sigma_B"*string("[cts/keV/kg/yr]")
        else

            res=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))
            bias=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))
            aff=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))

            for (idx,part) in enumerate(partitions)

                if (part_event_index[idx]!=0)
                    i_new = part_event_index[idx]
                    res[i_new]=Truncated(Normal(part.fwhm/2.355,part.fwhm_sigma/2.355),0,Inf)
                    bias[i_new] =Truncated(Normal(part.bias,part.bias_sigma),-Inf,Inf)
                    eff[i_new] =Truncated(Normal(part.eff_tot,part.eff_tot_sigma),0,1)
                    long_name = string(part.experiment)*" "*string(part.part_name)*" "*part.detector
                    append!(pretty_names[:σ],["Energy Resolution "*L"(\sigma)"*" "*long_name*" [keV]"])
                    append!(pretty_names[:𝛥],["Energy Scale Bias "*L"(\Delta)"*" - "*long_name*" [keV]"])
                    append!(pretty_names[:ε],["Efficiency "*L"(\varepsilon)"*" - "*long_name*""])
                end
            end

            # make some nice names for plotting
            for key in keys(distrB_multi)
                pretty_names[key]=string(key)*" [cts/keV/kg/yr]"
            end

            dis_B = distprod
            if bkg_only==false
                hd = BAT.HierarchicalDistribution(
                        v -> begin 
                        dict = (; (key =>LogNormal(log(v.B)-0.5*v.σB*v.σB,v.σB) for key in keys(distrB_multi))...)
                        BAT.NamedTupleDist(;dict...)
                        end,
                        BAT.NamedTupleDist(S=distrS,B=distrB,σB=0..1
                        , ε=eff, σ=res, 𝛥=bias
                        )
                ) 
            else
                hd = BAT.HierarchicalDistribution(
                        v -> begin 
                        dict = (; (key =>LogNormal(log(v.B)-0.5*v.σB*v.σB,v.σB) for key in keys(distrB_multi))...)
                        BAT.NamedTupleDist(;dict...)
                        end,
                        BAT.NamedTupleDist(B=distrB,σB=0..1
                        , ε=eff, σ=res, 𝛥=bias
                        )
                ) 
            end
            pretty_names[:B]="B [cts/keV/kg/yr]"
            pretty_names[:σB]=L"\sigma_B"*string("[cts/keV/kg/yr]")
            
        end

        return hd,pretty_names


    else 

        for key in keys(distrB_multi)
            pretty_names[key]=string(key)*" [cts/keV/kg/yr]"
        end

        dis_B = distprod
        hd = BAT.HierarchicalDistribution(
                v -> begin 
                dict = (; (key =>LogNormal(log(v.B)-0.5*v.σB*v.σB,v.σB) for key in keys(distrB_multi))...)
                BAT.NamedTupleDist(;dict...)
                end,
                BAT.NamedTupleDist(S=distrS,B=distrB,σB=0..1
                )
        )          
        pretty_names[:B]="B [cts/keV/kg/yr]"
        pretty_names[:σB]=L"\sigma_B"*string("[cts/keV/kg/yr]")


        return hd,pretty_names

    end

end

