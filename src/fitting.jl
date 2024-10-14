using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets
using TypedTables
using Plots
using Cuba

include("likelihood.jl")


function get_bkg_info(config)
    bkg_shape=:uniform
    bkg_shape_pars=nothing

    if (haskey(config["bkg"],"shape"))
        bkg_shape = Symbol(config["bkg"]["shape"]["name"])
        if (haskey(config["bkg"]["shape"],"pars"))
            bkg_shape_pars=config["bkg"]["shape"]["pars"]
        end
    end
    return bkg_shape,bkg_shape_pars
end


function norm_linear(x::Float64,p::NamedTuple,b_name::Symbol,fit_range)
    """
    Normalised linear function defined by (1+slope*(x-center)/260)/norm.
    Parameters
    ----------
        - slope::Real, the slope of the background
        - x::Real,     the x value to evaluate at
    """
        center = 1930
        range_l = [arr[1] for arr in fit_range]
        range_h = [arr[2] for arr in fit_range]
       
        sum_range = sum(range_h .- range_l)
        sum_range_sq =sum(range_h .^ 2 .- range_l .^ 2)
        slope = p[Symbol(string(b_name)*"_slope")]
    
        delta = range_h[end] - range_l[1]
        norm = sum_range * (1 - slope * center / delta) + slope * sum_range_sq / (2 * delta)
     
        return (1+slope*(x-center)/delta)/norm
    end



function norm_uniform(x::Real,p::NamedTuple,b_name::Symbol,fit_range)
    """
    Normalised linear function defined by (1+slope*(x-center)/260)/norm.
    Parameters
    ----------
        - x::Real,     the x value to evaluate at
    """    
    center = 1930
    range_l = [arr[1] for arr in fit_range]
    range_h = [arr[2] for arr in fit_range]

    norm =sum(range_h .- range_l)
    return 1/norm
 

end
 
function exp_stable(x::Float64)
    if (abs(x)<1E-6)
        return 1 + x + x^2/2+x^3/6
    else
        return exp(x)
    end
end
function norm_exponential(x::Float64,p::NamedTuple,b_name::Symbol,fit_range)
        """
        Normalised linear function defined by (1+slope*(x-center)/fit_range)/norm.
        Parameters
        ----------
            - slope::Real, the slope of the background
            - x::Real,     the x value to evaluate at
        """
        center = 1930
        range_l = [arr[1] for arr in fit_range]
        range_h = [arr[2] for arr in fit_range]

        centers=[center,center,center]
        # could be made faster?
        R = p[Symbol(string(b_name)*"_slope")]
        delta = range_h[end] - range_l[1]
        Rt=R/delta
        if (abs(Rt)>1E-6)
            norm = (-sum(exp_stable.(-Rt*(centers-range_l)))+sum(exp_stable.(-Rt*(centers-range_h))))/Rt
        else
            norm =sum(range_h .- range_l)
        end
       
        return exp_stable((x-center)*Rt)/norm
    
end




##############################################
##############################################
##############################################
function get_stat_blocks(partitions,events::Array{Vector{Float64}},part_event_index,fit_ranges;config,bkg_only)
"""
Function to retrieve useful pieces (prior, likelihood, posterior), also in saving values
"""
    settings=get_settings(config)

    
    if (haskey(config["bkg"],"correlated")) & (config["bkg"]["correlated"]["mode"]!="none")
        corr= true
        hier_mode =config["bkg"]["correlated"]["mode"]
        hier_range =config["bkg"]["correlated"]["range"]
    else
        corr=false
        hier_mode=nothing
        hier_range=nothing
    end

    bkg_shape,bkg_shape_pars = get_bkg_info(config)       
    
    prior,par_names=build_prior(partitions,part_event_index,config,settings,hierachical=corr,hierachical_mode=hier_mode,hierachical_range=hier_range,
                            bkg_shape=bkg_shape,shape_pars=bkg_shape_pars)
    @info "using a ",bkg_shape," bkg with ",bkg_shape_pars," parameters"
    @info "built prior"
    
    sqrt_prior=false
    s_max=nothing
    if bkg_only==false
        if (config["signal"]["prior"]=="sqrt")
            sqrt_prior=true
            s_max = Float64(config["signal"]["upper_bound"])
        end
    end
   
    likelihood = build_likelihood_looping_partitions(partitions, events, part_event_index,settings,sqrt_prior,s_max,fit_ranges,bkg_shape=bkg_shape)
    @info "built likelihood"
    
    posterior = PosteriorMeasure(likelihood, prior) 
    @info "got posterior"
    return prior,likelihood,posterior,par_names
end


function run_fit_over_partitions(partitions,events::Array{Vector{Float64}},part_event_index::Vector{Int}, config,fit_ranges)
"""
Function to run the fit looping over partitions
"""
    bkg_only = config["bkg_only"]
    prior,likelihood,posterior,par_names = get_stat_blocks(partitions,events,part_event_index,fit_ranges,config=config,bkg_only=bkg_only)

    Ns = Int(config["bat_fit"]["nsteps"])
    Nc = Int(config["bat_fit"]["nchains"])
    return bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = Ns, nchains = Nc)).result,prior,par_names
end
    

##############################################
##############################################
##############################################
function get_evidence(data,func,prior,method)
    likelihood =build_simple_likelihood(data,func)
    posterior = PosteriorMeasure(likelihood, prior)
    
    return bat_integrate(posterior,method)
end






##############################################
##############################################
##############################################
function get_qbb_posterior(fit_function,samples)
    qbb=[]
    
    for samp in samples
        v=samp.v
        weight=samp.weight
        for w in 1:1:weight
            append!(qbb,fit_function(NamedTuple(v),2039.0))
        end
    end
    return qbb
end

##############################################
##############################################
##############################################
function get_par_posterior(samples,par;idx)
    
    pars=[]
    
    for samp in samples
        v=samp.v
        weight=samp.weight
        for w in 1:1:weight
            if (idx==nothing)
                append!(pars,v[par])
            else
                append!(pars,v[par][idx])
            end
        end
    end
   
    return pars
end


##############################################
##############################################
##############################################
function get_n_posterior(samples)
    ns=[]
    
    for samp in samples
        v=samp.v
        weight=samp.weight
        for w in 1:1:weight
            append!(ns,v.n)
        end
    end
    return ns
end




struct LogFlat <: ContinuousUnivariateDistribution
        a::Float64
        b::Float64
    end

Distributions.support(::LogFlat) = (0.0, Inf)
Distributions.pdf(d::LogFlat, x) = (x > d.a) && (x<d.b) ? 1/(x*(log(d.b)-log(d.a))) : 0
Distributions.rand(d::LogFlat) = exp(rand()*(log(d.b)-log(d.a))+log(d.a))
Distributions.logpdf(d::LogFlat, x::Float64) = x > d.a && x<d.b ? -log(x)-log(log(d.b)-log(d.a)) : -Inf
Distributions.cdf(d::LogFlat, x::Float64) =  x < d.a ? 0 : x<d.b ? (log(x)-log(d.a))/(log(d.b)-log(d.a)) : 1
Distributions.params(d::LogFlat) = (d.a,d.b)

    
