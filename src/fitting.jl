using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets
using TypedTables
using Plots
using Cuba

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
function build_likelihood(data,func)
    likelihood = let d=data,f=func

        
        total_counts= length(d)
        logfuncdensity(function (params)
                norm =float(params.b)
                ll_value=0
                ll_value = logpdf(Poisson(norm+eps(norm)), total_counts)
                for energy in d
                    
                    ll_value+=log((f(params,energy)+eps(norm))/(eps(norm)+norm))
                end

         
            return ll_value
        end)
    end
end



##############################################
##############################################
##############################################
function build_likelihood_signal_background(data,func_bkg,func_sig,fwhm,is_uniform = true)
    """
    Build the likelihood as a sum of normalised functions
    """
    likelihood = let d=data,fs=func_sig,fb=func_bkg

        # convert fwhm to sigma
        sigma = fwhm/2.355
        total_counts= length(d)

        # define the likelihood
        logfuncdensity(function (params)

                # normalisation
                norm =float(params.b)+float(params.n)
                
                ll_value=0

                # poisson term
                ll_value = logpdf(Poisson(norm+eps(norm)), total_counts)

                # extract slope (either 0 or from parameters)
                if (is_uniform)
                    slope = 0
                else
                    slope = params.s
                end

                # loop over energies
                for energy in d

                    #could be edited for systematics on sigma / mu
                    model_s = float(params.n)*fs(sigma,2039,energy)
                    model_b = float(params.b)*fb(slope,energy)
                    model = float(model_s)+float(model_b)+eps(model_s)+eps(model_b)
                    
                    norm = norm+eps(norm)
                    
                    ll_value+=log((model_s+model_b)/norm)
                    
                end

         
            return ll_value
        end)
    end
end



##############################################
##############################################
##############################################
function run_fit(data,func,prior)
    likelihood = build_likelihood(data,func)
    posterior = PosteriorMeasure(likelihood, prior)

    return bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)).result
end



##############################################
##############################################
##############################################
function get_evidence(data,func,prior,method)
    likelihood =build_likelihood(data,func)
    posterior = PosteriorMeasure(likelihood, prior)
    
    return bat_integrate(posterior,method)
end



##############################################
##############################################
##############################################
function run_fit_signal(data,prior,is_uniform=true)
    likelihood =build_likelihood_signal_background(data,norm_uniform,norm_gauss,3,is_uniform)
    print(likelihood)
    print(prior)
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



