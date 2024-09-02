using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using BAT, DensityInterface, IntervalSets
using TypedTables
using Plots
using Cuba

# define some constants
Qbb = 2039.06 # keV
N_A = 6.022E23
m_76 = 75.6E-3 # kg/mol
deltaE = 240 # keV



# the function is called in a loop of partitions (and experiments)
function build_likelihood_per_partition(idx_k, part_k::Table, events, p)

    # free parameters: signal (S), background (B), energy bias (biask) and resolution per partition (resk)
    ll_value = 0

    model_s_k = log(2) * N_A * part_k.livetime * part_k.mass * part_k.eff_k.val * p.S / m_76
    model_b_k = deltaE * part_k.livetime * part_k.mass * p.B
    model_tot_k = model_b_k + model_s_k

    # loop over events in the analysis window (we already checked for presence of events for a given partition k)
    for i in part_k.events
        """
        # old style
        for i in events

            # for a fixed partition, load all energy events falling in that partition (void if no match is found)
            events_in_k = events.i[
                (events.timestamp .>= part_k.span_in_utc_s[1]) .&&
                (events.timestamp .<= part_k.span_in_utc_s[2]) .&&
                (events.detector == part_k.detector)
            ]
        """

        # poisson term
        ll_value += logpdf(Poisson(model_tot_k), length(events_in_k)) # + alpha term ???

        for evt_energy in events_in_k

            term1 = model_b_k / deltaE # background
            term2 = model_s_k * pdf(Normal(Qbb + p.bias[idx_k], p.res[idx_k]), evt_energy) # signal

            ll_value += log( (term1 + term2) / model_tot_k ) 
        end

        """ # bring outside
        # retrieve val/unc for adding pull terms
        res_k = part_k.fwhm_in_keV.val / 2.355 # resolution is saved as fwhm for l200
        sigma_res_k = part_k.fwhm_in_keV.unc / 2.355
        bias_k = part_k.energy_bias_in_keV.val
        sigma_bias_k = part_k.energy_bias_in_keV.unc

        term3 = pdf(Normal(bias_k, sigma_bias_k), p.bias[idx_k]) # pull term for energy bias
        term4 = pdf(Normal(res_k, sigma_res_k), p.res[idx_k]) # pull term for energy resolution
        ### TO DO - alpha efficiency term
        ll_value += log( (term3 * term4) / model_tot_k ) 
        """

    end

    return ll_value
end

# Tuple{Real, Real, Vector{Real}, Vector{Real}}
ModelParameters = NamedTuple{(:S, :B, :bias, :res)}
function build_likelihood_looping_partitions(partitions, events)
    DensityInterface.logfuncdensity(
        p::ModelParameters -> begin
            total_ll = 0.0

            for (idx_k, part_k) in enumerate(partitions)
                total_ll += build_likelihood_per_partition(idx_k, part_k, events, p)
            end

            return total_ll
        end
    )
end

##############################################
##############################################
##############################################
function build_simple_likelihood(data,func)
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
function build_simple_likelihood_signal_background(data,func_bkg,func_sig,fwhm,is_uniform = true)
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
                    model_s = float(params.n)*fs(sigma,Qbb,energy)
                    model_b = float(params.b)*fb(slope,energy)
                    model = float(model_s)+float(model_b)+eps(model_s)+eps(model_b)
                    
                    norm = norm+eps(norm)
                    
                    ll_value+=log((model_s+model_b)/norm)
                    
                end

         
            return ll_value
        end)
    end
end



