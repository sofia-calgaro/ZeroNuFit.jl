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


function build_likelihood_zero_obs_evts(part_k, p;stat_only=false)

    ll_value = 0
    
    model_s_k = log(2) * N_A * part_k.exposure * part_k.eff_tot * p.S / m_76
    model_b_k = deltaE * part_k.exposure * p.B
    model_tot_k = model_b_k + model_s_k

    ll_value += -(model_tot_k+eps(model_tot_k)) # + alpha term ???
    
    return ll_value
end

function build_likelihood_per_partition(idx_k, part_k, events_k, p;stat_only=false)

    # free parameters: signal (S), background (B), energy bias (biask) and resolution per partition (resk)
    ll_value = 0
    
    model_s_k = log(2) * N_A * part_k.exposure * part_k.eff_tot * p.S / m_76
    model_b_k = deltaE * part_k.exposure * p.B
    model_tot_k = model_b_k + model_s_k

    ll_value += logpdf(Poisson(model_tot_k+eps(model_tot_k)), length(events_k)) # + alpha term ???
    
    for i in events_k
        for evt_energy in events_k

            term1 = model_b_k / deltaE # background

            if (stat_only==true)
                term2 = model_s_k * pdf(Normal(Qbb + part_k.bias, part_k.fwhm/2.355), evt_energy) # signal (fixed nuisance)
            else
                term2 = model_s_k * pdf(Normal(Qbb + p.bias[idx_k], p.res[idx_k]), evt_energy) # signal (free nuisance)
            end
            ll_value += log( (term1 + term2)+eps(term1+term2)) - log(model_tot_k+eps(model_tot_k)) 
        end

    end
    
    return ll_value
end

# Tuple{Real, Real, Vector{Real}, Vector{Real}}
ModelParameters = NamedTuple{(:S, :B, :bias, :res)}
function build_likelihood_looping_partitions(partitions, events;stat_only=false)
    return DensityInterface.logfuncdensity( function(p)
            total_ll = 0.0

            for (idx_k, part_k) in enumerate(partitions)
                
                if events[idx_k] != Any[]
                    total_ll += build_likelihood_per_partition(idx_k, part_k, events[idx_k], p, stat_only=stat_only)
                else
                    # no events are there for a given partition
                    total_ll += build_likelihood_zero_obs_evts(part_k, p)
                end
            end
            return total_ll
        end
    )
end

function build_prior(partitions;config,stat_only=false)
    res=[]
    bias=[]
    for part in partitions
        append!(res,[Normal(part.fwhm/2.355,part.fwhm_sigma/2.355)])
        append!(bias,[Normal(part.bias,part.bias_sigma)])
    end
    if (stat_only==false)
        return distprod(S=0..config["upper_signal"],B=0..config["upper_bkg"], res= fill(0..10, length(partitions)),
    bias= fill(-2..2,length(partitions)))
    else 
        distprod(S=0..config["upper_signal"],B=0..config["upper_bkg"])
    end
    
    #res=distprod(res),bias=distprod(bias))


end

