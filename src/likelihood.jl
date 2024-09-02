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
function build_likelihood_per_partition(idx_k, part_k::Table, events_k, p)

    # free parameters: signal (S), background (B), energy bias (biask) and resolution per partition (resk)
    ll_value = 0

    model_s_k = log(2) * N_A * part_k.exposure * part_k.eff_tot * p.S / m_76
    model_b_k = deltaE * part_k.exposure * p.B*part_k.eff_tot
    model_tot_k = model_b_k + model_s_k

    # loop over events in the analysis window (we already checked for presence of events for a given partition k)
    
        
    ll_value += logpdf(Poisson(model_tot_k), length(events_k)) # + alpha term ???
    for i in events_k
        for evt_energy in events_in_k

            term1 = model_b_k / deltaE # background
            term2 = model_s_k * pdf(Normal(Qbb + p.bias[idx_k], p.res[idx_k]), evt_energy) # signal

            ll_value += log( (term1 + term2) / model_tot_k ) 
        end

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
                total_ll += build_likelihood_per_partition(idx_k, part_k, events[idx_k], p)
            end

            return total_ll
        end
    )
end

