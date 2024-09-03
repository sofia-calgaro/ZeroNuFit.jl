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
   
    norm = sum(range_h .- range_l) * (1 - slope * center / 260) + slope * sum(range_h .^ 2 .- range_l .^ 2) / (2 * 260)
    (1+slope*(x-center)/260)/norm
end


##############################################
##############################################
##############################################
function norm_gauss(sigma::Real,mu::Real,x::Real)
    pdf(Normal(mu,sigma), x)
end


##############################################
##############################################
##############################################
function fit_function_linear(p::NamedTuple{(:b, :s)}, x::Real)
    p.b*norm_uniform(p.s,x)
end
function fit_function_uniform(p::NamedTuple{(:b,)},x::Real)
    p.b*norm_uniform(0,x)
end



##############################################
##############################################
##############################################
function run_fit_over_partitions(partitions,events;func,config,stat_only)

    prior=build_prior(partitions,config=config,stat_only=stat_only)
    @info "build prior"
    likelihood = build_likelihood_looping_partitions(partitions, events,stat_only=stat_only)
    posterior = PosteriorMeasure(likelihood, prior) 
    @info "got posterior"

    return bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4),).result
end



##############################################
##############################################
##############################################
function run_fit(data,func,prior)
    likelihood = build_simple_likelihood(data,func)
    posterior = PosteriorMeasure(likelihood, prior)

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
function run_fit_signal(data,prior,is_uniform=true)
    likelihood =build_simple_likelihood_signal_background(data,norm_uniform,norm_gauss,3,is_uniform)
   
    posterior = PosteriorMeasure(likelihood, prior)

    return bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)).result
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



