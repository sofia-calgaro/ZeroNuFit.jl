using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets
using TypedTables
using Plots
using Cuba

include("likelihood.jl")

# we need to define some building blocks
#TODO: move to config
center = 1930
range_l = [1930, 2109, 2124]
range_h = [2099, 2114, 2190]

##############################################
##############################################
##############################################
function norm_uniform(slope::Real,x::Real)
"""
Normalised linear function defined by (1+slope*(x-center)/260)/norm.
Parameters
----------
    - slope::Real, the slope of the background
    - x::Real,     the x value to evaluate at
"""
    norm = sum(range_h .- range_l) * (1 - slope * center / 260) + slope * sum(range_h .^ 2 .- range_l .^ 2) / (2 * 260)
    (1+slope*(x-center)/260)/norm
end


##############################################
##############################################
##############################################
function norm_gauss(sigma::Real,mu::Real,x::Real)
"""
Normalised gaussian function
Parameters
----------
    - sigma::Real (the std of the Gaussian)
    - mu::Real    (the mean)
    - x::Real     (the x value to evaluate at)
"""
    pdf(Normal(mu,sigma), x)
end





##############################################
##############################################
##############################################
function get_stat_blocks(partitions,events::Array{Vector{Float64}},part_event_index;config,bkg_only)
"""
Function to retrieve useful pieces (prior, likelihood, posterior), also in saving values
"""
    settings=get_settings(config)

    # check if the key 'correlated' exists
    if !haskey(config["bkg"], "correlated")
        throw(ArgumentError("Key 'correlated' not found for the background parameter in the configuration JSON! Exit here"))
    end
    
    corr= config["bkg"]["correlated"]
    
    prior,par_names=build_prior(partitions,part_event_index,config,settings,hierachical=corr)
    
    @info "built prior"
    
    sqrt_prior=false
    s_max=nothing
    if bkg_only==false
        if (config["signal"]["prior"]=="sqrt")
            sqrt_prior=true
            s_max = config["signal"]["upper_bound"]
        end
    end
    likelihood = build_likelihood_looping_partitions(partitions, events, part_event_index,settings,sqrt_prior)
    @info "built likelihood"
    
    posterior = PosteriorMeasure(likelihood, prior) 
    @info "got posterior"
    return prior,likelihood,posterior,par_names
end


function run_fit_over_partitions(partitions,events::Array{Vector{Float64}},part_event_index::Vector{Int}, config)
"""
Function to run the fit looping over partitions
"""
    bkg_only = config["bkg_only"]
    prior,likelihood,posterior,par_names = get_stat_blocks(partitions,events,part_event_index,config=config,bkg_only=bkg_only)

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

    
