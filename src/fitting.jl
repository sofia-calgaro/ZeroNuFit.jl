using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets
using TypedTables
using Plots
using Cuba

include("likelihood.jl")

# we need to define some building blocks
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
function run_fit_over_partitions(partitions,events,part_event_index;config,stat_only)
"""
FUnction to run the fit looping over partitions
"""
    println(part_event_index)
    prior=build_prior(partitions,part_event_index,config=config,stat_only=stat_only)
    @info "build prior"
    likelihood = build_likelihood_looping_partitions(partitions, events, part_event_index,stat_only=stat_only)
    posterior = PosteriorMeasure(likelihood, prior) 
    @info "got posterior"

    return bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)).result
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



