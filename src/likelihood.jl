using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets

using TypedTables
using Plots,LaTeXStrings
using Cuba
using OrderedCollections



function get_bkg_pdf(bkg_shape::Symbol,evt_energy::Float64,p::NamedTuple,b_name::Symbol,fit_range)
    if (bkg_shape==:uniform)
        return norm_uniform(evt_energy,p,b_name,fit_range)
    elseif (bkg_shape==:linear)
        return norm_linear(evt_energy,p,b_name,fit_range)
    elseif (bkg_shape==:exponential)
        return norm_exponential(evt_energy,p,b_name,fit_range)
    else
        @error "bkg shape",bkg_shape," is not yet implememnted"
        exit(-1)
    end

end
function get_energy_scale_pars(part_k::NamedTuple,p::NamedTuple,settings::Dict,idx_part_with_events)
    """ 
    Get the resolution and bias
    """
    if (settings[:energy_scale_fixed]==true || idx_part_with_events==0)
        reso = part_k.fwhm/2.355
        bias =part_k.bias

    elseif (settings[:energy_scale_correlated]==true)
        energy_reso_group = part_k.energy_reso_name
        energy_bias_group = part_k.energy_bias_name
        reso = part_k.fwhm/2.355+p[energy_reso_group]*part_k.fwhm_sigma/2.355
        bias = part_k.bias+p[energy_bias_group]*part_k.bias_sigma
        
    else
        reso = p.σ[idx_part_with_events]
        bias =p.𝛥[idx_part_with_events]
    end

return reso,bias

end

function get_mu_s_b(p::NamedTuple,part_k::NamedTuple,idx_part_with_events::Int,settings::Dict,fit_range)
    """
    Get the expected number of signal and background counts in a partition
    """
    N_A = 6.022E23
    m_76 = 75.92E-3 # kg/mol
    sig_units =1e-27 # signal is in units of this
    
    deltaE = sum([arr[2]-arr[1] for arr in fit_range])
    eff= nothing
    
    # logic for efficiency it can be either correlated, uncorrelated or fixed
    if settings[:bkg_only]==true 
        eff =0

    elseif (settings[:eff_correlated] == true)
        eff_group = part_k.eff_name
        eff = part_k.eff_tot + p[eff_group] * part_k.eff_tot_sigma

    elseif (idx_part_with_events!=0 && 
            settings[:eff_correlated]==false &&
            settings[:eff_fixed]==false)

        eff =p.ε[idx_part_with_events] 
    else 
        eff =part_k.eff_tot
    end

    if (settings[:bkg_only]==false)
        model_s_k = log(2) * N_A * part_k.exposure * (eff) * (p.S*sig_units) / m_76
    else
        model_s_k=0
    end
    
    b_name = part_k.bkg_name
    model_b_k = deltaE * part_k.exposure * p[b_name]

    return model_s_k,model_b_k
end


function build_likelihood_zero_obs_evts(part_k::NamedTuple, p::NamedTuple,settings::Dict,fit_range)
    """
    Function to calculate the partial likelihood for a partition with 0 events
    """

    ll_value = 0
    model_s_k,model_b_k = get_mu_s_b(p,part_k,0,settings,fit_range)
    model_tot_k = model_b_k + model_s_k

    ll_value += -(model_tot_k+eps(model_tot_k)) 
    
    return ll_value
end


function build_likelihood_per_partition(idx_k::Int, idx_part_with_events::Int,part_k::NamedTuple, events_k::Vector{Float64}, 
    p::NamedTuple,settings::Dict,bkg_shape::Symbol,fit_range)
"""
Function which computes the partial likelihood for a single data partiton
free parameters: signal (S), background (B), energy bias (biask) and resolution per partition (resk)
"""
    Qbb = 2039.06 # keV

    ll_value = 0

    model_s_k,model_b_k =   get_mu_s_b(p,part_k,idx_part_with_events,settings,fit_range)
   
    model_tot_k = model_b_k + model_s_k
    
    # constrain λ not to be negative
    if model_tot_k < 0 
        λ = eps(0.)
    else
        λ = model_tot_k+eps(model_tot_k)
    end

    ll_value += logpdf(Poisson(λ), length(events_k))

    for evt_energy in events_k
      
        term1 = model_b_k * get_bkg_pdf(bkg_shape,evt_energy,p, part_k.bkg_name ,fit_range)

        if (settings[:bkg_only]==false)

            # get the correct reso and bias (
            reso,bias = get_energy_scale_pars(part_k,p,settings,idx_part_with_events)
            term2 = model_s_k * pdf(Normal(Qbb - bias, reso), evt_energy) 
        else
            term2 =0
        end
        
        ll_value += log( (term1 + term2)+eps(term1+term2)) - log(model_tot_k+eps(model_tot_k)) 
    
    end
    
    return ll_value
end



function build_likelihood_looping_partitions(partitions::TypedTables.Table, events::Array{Vector{Float64}},
    part_event_index::Vector{Int},settings::Dict,sqrt_prior::Bool,s_max::Union{Float64,Nothing},fit_ranges;bkg_shape::Symbol=:uniform)
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
                    total_ll += build_likelihood_per_partition(idx_k,part_event_index[idx_k], part_k, events[idx_k], p, settings,bkg_shape,fit_ranges[part_k.fit_group])
                else
                    # no events are there for a given partition
                    total_ll += build_likelihood_zero_obs_evts(part_k, p,settings,fit_ranges[part_k.fit_group])
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




function generate_data(samples::BAT.DensitySampleVector,partitions::TypedTables.Table,part_event_index::Vector{Int},settings::Dict,fit_ranges;
    best_fit::Bool=false,seed=nothing,bkg_only=false)
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
        model_s_k,model_b_k = get_mu_s_b(p,part_k,idx_part_with_events,settings,fit_ranges[part_k.fit_group])

        n_s = rand(Poisson(model_s_k))
        n_b = rand(Poisson(model_b_k))
        events =generate_disjoint_uniform_samples(n_b)
        if (bkg_only == false)
            for i in 1:n_s

                reso,bias = get_energy_scale_pars(part_k,p,settings,idx_part_with_events)

                append!(events,rand(Normal(Qbb - bias, reso)))
                
        
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

function build_prior(partitions,part_event_index,config,settings::Dict;hierachical=false,hierachical_mode=nothing,hierachical_range=nothing,bkg_shape=:uniform,shape_pars=nothing)
"""
Builds the priors for use in the fit
----------
Parameters
    - partitions:Table of the partition info
    - config: the Dict of the fit config
    - nuis_prior; true if we want to include priors for nuisance parameters (bias, res, eff)
"""

    # bkg indexs
    list_names = partitions.bkg_name
    unique_list=unique(list_names)
    bkg_names=[Symbol(name) for name in unique_list]
    
   

    distrS, distrB = get_signal_bkg_priors(config)
    distrB_multi=OrderedDict(Symbol(bkg_name)=>distrB for bkg_name in bkg_names)

    pretty_names =OrderedDict(:S=>string("S [")*L"10^{-27}"*string("yr")*L"^{-1}"*string("]"),
        :α=>L"\alpha_{\varepsilon}",
        :αr=>L"\alpha_{r}",
        :αb=>L"\alpha_{b}",
        :ε=>[],
        :σ=>[],
        :𝛥=>[])

    for key in keys(distrB_multi)
        pretty_names[key]=string(key)*" [cts/keV/kg/yr]"
    end
    

    # create priors one by one
    ### SIGNAL PRIOR
    
    priors=OrderedDict()
    if (settings[:bkg_only]==false)
        priors[:S] =distrS
        @info "entered to add S prior"
    end

    ### EFF prior

    if (settings[:eff_fixed]==false  && settings[:eff_correlated]==true)

        all_eff_tot = partitions.eff_tot
        all_eff_tot_sigma = partitions.eff_tot_sigma
        ratio = - all_eff_tot ./ all_eff_tot_sigma 
        α_min = maximum(ratio)

        list_names = partitions.eff_name
        unique_list=unique(list_names)
        for name in unique_list
            priors[Symbol(name)]=Truncated(Normal(0,1),α_min,Inf)
            pretty_names[Symbol(name)]=L"\alpha_{\varepsilon} ("*split(String(name),"_")[2]*")"
        end

      
    else
        eff=Vector{Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}}(undef,maximum(part_event_index))
    

        for (idx,part) in enumerate(partitions)

            if (part_event_index[idx]!=0)
                i_new = part_event_index[idx]
                
                eff[i_new] =Truncated(Normal(part.eff_tot,part.eff_tot_sigma),0,1)
                long_name = string(part.experiment)*" "*string(part.part_name)*" "*part.detector
                append!(pretty_names[:ε],["Efficiency "*L"(\varepsilon)"*" - "*long_name*""])
            end
        end
        priors[:ε]=eff

    end

    ### ENERGY scale prior

    if (settings[:energy_scale_fixed]==false && settings[:energy_scale_correlated]==true)
        all_fwhm = partitions.fwhm
        all_fwhm_sigma = partitions.fwhm_sigma
        ratio = - all_fwhm ./ all_fwhm_sigma 
        αr_min = maximum(ratio)

        list_names = partitions.energy_reso_name
        unique_list=unique(list_names)
        for name in unique_list
            priors[Symbol(name)]=Truncated(Normal(0,1),αr_min,Inf)
            pretty_names[Symbol(name)]=L"\alpha_{r} ("*split(String(name),"_")[2]*")"

        end


        list_names = partitions.energy_bias_name
        unique_list=unique(list_names)
        for name in unique_list
            priors[Symbol(name)]=Truncated(Normal(0,1),-Inf,Inf)
            pretty_names[Symbol(name)]=L"\alpha_{b} ("*split(String(name),"_")[2]*")"

        end

    else
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
        priors[:σ]=res
        priors[:𝛥]=bias

    end

    ## bkg shape priors
    if shape_pars!=nothing
        
        for par in keys(shape_pars)
            name = par
            prior = shape_pars[par]
          
            for bkg_name in bkg_names
                priors[Symbol(string(bkg_name)*"_"*name)]=Uniform(prior[1],prior[2])
                pretty_names[Symbol(string(bkg_name)*"_"*name)]=string(bkg_name)*"_"*name
            end

        end
    
    end
    
    ## BKG prior
    if (hierachical==false)
        for (key,item) in distrB_multi
            priors[key]=item
        end

        return distprod(;priors...),pretty_names
    else

        if (hierachical_mode=="lognormal")
            hd = BAT.HierarchicalDistribution(
                v -> begin 
                dict = (; (key =>LogNormal(log(v.B)-0.5*v.σB*v.σB,v.σB) for key in keys(distrB_multi))...)
                return distprod(;dict...,priors...)
                end,
                distprod(B=0..1,σB=Uniform(hierachical_range[1],hierachical_range[2]))
                
                )
        elseif (hierachical_mode=="normal")
            hd = BAT.HierarchicalDistribution(
                v -> begin 
                dict = (; (key =>Normal(log(v.B)-0.5*v.σB*v.σB,v.σB) for key in keys(distrB_multi))...)
                return distprod(;dict...,priors...)
                end,
                distprod(B=0..1,σB=Uniform(hierachical_range[1],hierachical_range[2]))
                
                )
        else
            @error "hierachical (correlated) bkg mode $hierachical_mode is not know"
            exit(-1)
        end

        pretty_names[:B]="B [cts/keV/kg/yr]"
        pretty_names[:σB]=L"\sigma_B"*string("[cts/keV/kg/yr]")

        return hd,pretty_names
    end
end

function print_names(priors,pretty_names)

    for (key, value) in pairs(priors)
        @info "added $key  to the model"
        if (pretty_names[key] isa Vector)
            for (idx,pn) in enumerate(pretty_names[key])
                 @info "added $key[$idx] ($pn)  to the model"
            end
        else
            pn = pretty_names[key]
            @info "added $key ($pn) to the model"

        end

    end

end

