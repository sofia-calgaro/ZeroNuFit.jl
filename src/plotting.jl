using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using Plots
using BAT, DensityInterface, IntervalSets
using TypedTables
using Cuba
using Base.Filesystem
using PDFmerger: append_pdf!

default(
    framestyle=:box,               # Grid line transparency
    background_color = :white   ,       # Background color of the plot,
    titlefontsize=12,     # Global title font size
    guidefontsize=12,     # Global axis label font size
    tickfontsize=12,      # Global tick label font size
    legendfontsize=8     # Global legend font size
)
# define some constants -TODO add to config
Qbb = 2039.06 # keV
N_A = 6.022E23
m_76 = 75.6E-3 # kg/mol
deltaE = 240 # keV
sig_units =1e-27 # signal is in units of this

##############################################
##############################################
##############################################
function constant(x)
    if (x<2124 && x>2114)
        0
    elseif (x<2109 && x>2099)
        0
    else   
        1/240 
    end
end


function fit_model(p, part_k, x, is_in_part, part_idx)
    
    model_s_k = log(2) * N_A * part_k.exposure * (part_k.eff_tot + p.α * part_k.eff_tot_sigma) * (p.S*sig_units) / m_76
    model_b_k = deltaE * part_k.exposure * p.B

    term1 = model_b_k / deltaE
    
    if is_in_part
        term2 = model_s_k * pdf(Normal(Qbb + p.𝛥[part_idx], p.σ[part_idx]), x)
    else
        term2 = model_s_k * pdf(Normal(Qbb + part_k.bias, part_k.fwhm/2.355), x)
    end
    
    return term1 + term2 
    
end


##############################################
##############################################
##############################################
function plot_data(hist::Histogram,name,partitions,part_event_index,pars,samples,plotflag)
"""
Function to plot events in the Qbb analysis window and BAT fit results
"""
    
    counts=sum(hist.weights)
    p = plot() 

    ymax = 1.5
    bin_edges = hist.edges[1]
    
    plot!(
       p, hist,
        st = :steps, label = "Data",
        title ="$name",
        xlabel="Energy (keV)",
        ylabel="Counts/ (1 keV)",
        ylim=(0,ymax),
        xlim=(1930,2190),
        color="dark blue",
        fill=true,
        framestyle = :box
        
    )
    
    #plot fit model
    function find_a_name(params,x)
        model = 0
        for (idx_k, part_k) in enumerate(partitions)
            if part_event_index[idx_k]!=0
                idx_k_with_events=part_event_index[idx_k]
                model += fit_model(params, part_k, x, true, part_event_index[idx_k])*diff(bin_edges)[1]
            else
                model += fit_model(params, part_k, x, false, nothing)*diff(bin_edges)[1]
            end
        end
        return model
    end
    
    if plotflag["bandfit_and_data"]
        plot!(p,1930:0.1:2190, find_a_name, samples, alpha=0.4,median=false,globalmode=false,fillalpha=0.3) #TO DO: take only some samples
        best_fit_pars = BAT.mode(samples)
        plot!(p,1930:0.1:2190,x -> find_a_name(best_fit_pars,x),label="Fit",lw=2,color="red")
    else
        plot!(p,1930:0.1:2190, find_a_name, samples, alpha=0.4,median=false,globalmode=true,fillalpha=0.3)
    end
    
    # exclude the gamma lines
    shape_x = [2114,2114,2124,2124]
    shape_x2 = [2099,2099,2109,2109]
    shape_y=[0,1.5,1.5,0]
    
    plot!(p, shape_x, shape_y, fillalpha=0.3,fill=true, line=false,fillcolor=:blue, label=false)
    plot!(p, shape_x2, shape_y, fillalpha=0.3,fill=true, line=false,fillcolor=:blue, label=false)
    
    return p
end



##############################################
##############################################
##############################################
function plot_fit_and_data(partitions, events, part_event_index, samples, pars, output, plotflag)
    
    # create histo with energies 
    energies = []
    for (idx_k, part_k) in enumerate(partitions)
        if events[idx_k] != Any[]
            for energy in events[idx_k]
                append!(energies, events[idx_k])
            end
        end
    end
    hist_data = append!(Histogram(1930:1:2190), energies)
    
    p_fit = plot_data(hist_data,"",partitions,part_event_index,pars,samples,plotflag)
    
    savefig(joinpath(output, "plots/fit_over_data.pdf"))
    
end



##############################################
##############################################
##############################################
function plot_marginal_distr(partitions,samples,pars,output;priors=nothing)    
"""
Function to plot 1D and 2D marginalized distributions (and priors)
"""
    
    name = split(output, "output/")[end]
    first_sample = samples.v[1]
    unshaped_samples, f_flatten = bat_transform(Vector, samples)
    @debug "Unshaped samples:", bat_report(unshaped_samples)
    
    # remove old file
    if isfile(joinpath(output, "plots/marg_posterior.pdf"))
        Filesystem.rm(joinpath(output, "plots/marg_posterior.pdf"),force=true)
    end
    
    # 1D posteriors
    ct = 1
    for par in pars
        par_entry = first_sample[par]
        
        if length(par_entry) == 1
            post = get_par_posterior(samples,par,idx=nothing)
            if (par==:S || par==:B)
                mini=0
            else
                mini=minimum(post)
            end

            p=plot(
            samples, par,
            mean = false, std = false, globalmode = true, marginalmode = true,
            nbins = 200, xlim=(mini,maximum(post))
            ) 
            x=range(mini, stop=maximum(post), length=1000)

            # plot prior
            if priors!=nothing
                y=pdf(priors[par],x)
                plot!(x,y,label="prior",color="grey")
            end

            # TO DO: add a way to constrain the posterior in [0; max from config] or [0; right-est entry on the x axis for signal]
            savefig(p,"temp.pdf")
            append_pdf!(joinpath(output, "plots/marg_posterior.pdf"), "temp.pdf", cleanup=true)
            ct += 1
            
        # multivariate parameters    
        else
            for idx in 1:length(par_entry) 
                post = get_par_posterior(samples,par,idx=idx)

                xlab = string("$(par)[$(idx)]")
                ylab = string("P($(par)[$(idx)])")
                
                p=plot(
                unshaped_samples, ct,
                mean = false, std = false, globalmode = true, marginalmode = true,
                nbins = 200, xlabel = xlab, ylabel = ylab, xlim=(minimum(post),maximum(post))
                )
                
                savefig(p,"temp.pdf")
                append_pdf!(joinpath(output, "plots/marg_posterior.pdf"), "temp.pdf", cleanup=true)
                ct += 1
            end
        end
    end
    
    # 2D posteriors
    plot(
        samples,
        mean = false, std = false, globalmode = false, marginalmode = true,
        nbins = 200
    )
    savefig(joinpath(output, "plots/all_marg_posterior_2D.pdf"))
        
end



##############################################
##############################################
##############################################
function plot_qbb_comp(qbb,name)
    p=nothing
    for (q,n) in zip(qbb,name)
        hist = append!(Histogram(0:10/240. /100. :30/240.),q)
        if (p==nothing)
            func=plot
        else
            func=plot!
        end
         p=func(
           hist,
            st = :steps, label = "$n",

            xlabel="b(Qbb) [cts/keV]",
            ylabel="prob [arb.]",
         
            fill=true,
            fillalpha=0.3,
            framestyle = :box)
              
    end
  
end



##############################################
##############################################
##############################################
function plot_n_comp(ns,name)
    p=nothing
    for (q,n) in zip(ns,name)
        hist = append!(Histogram(0:0.1:10),q)
        if (p==nothing)
            func=plot
        else
            func=plot!
        end
         p=func(
           hist,
            st = :steps, label = "$n",

            xlabel="n [cts]",
            ylabel="prob [arb.]",
         
            fill=true,
            fillalpha=0.3,
            framestyle = :box)
              
    end
  
end