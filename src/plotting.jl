using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using Plots
using BAT, DensityInterface, IntervalSets
using TypedTables
using Cuba

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



##############################################
##############################################
##############################################
function plot_data(hist::Histogram,name,func=nothing,samples=nothing,fitfunction=nothing)
    counts=sum(hist.weights)
    p = plot() 

    ymax = 1.5
    bin_edges = hist.edges[1]
    function binfitfunction(x,params)
        fitfunction(x,params)*diff(bin_edges)[1]
    end
    
    if (samples!=nothing)
        plot!(p,1930:0.1:2190, binfitfunction, samples,alpha=0.4,median=false,globalmode=false,fillalpha=0.3)
    end
    if (func !=nothing)
        plot!(p,1930:0.1:2190,x -> diff(bin_edges)[1]*counts*func(x),label="constant",lw=2,color="red")
    end
    plot!(
       p, hist,
        st = :steps, label = "Data",
        title ="$name data",
        xlabel="Energy [keV]",
        ylabel="counts/2 keV",
        ylim=(0,ymax),
        xlim=(1930,2190),
        color="dark blue",
        fill=true,
        framestyle = :box
        
    )
  
    shape_x = [2114,2114,2124,2124]
    shape_x2 = [2099,2099,2109,2109]

    shape_y=[0,1.5,1.5,0]
    
    # exclude the gamma lines
    plot!(p, shape_x, shape_y, fillalpha=0.3,fill=true, line=false,fillcolor=:blue, label=false)
    plot!(p, shape_x2, shape_y, fillalpha=0.3,fill=true, line=false,fillcolor=:blue, label=false)
    
    return p
end



##############################################
##############################################
##############################################

function make_plots(partitions,samples,pars,output)    
    
    name = split(output, "output/")[end]
    first_sample = samples.v[1]
    unshaped_samples, f_flatten = bat_transform(Vector, samples)
    @debug "Unshaped samples:", bat_report(unshaped_samples)
    
    # marginalized posterior for each parameter
    ct = 1
    for par in pars
        par_entry = first_sample[par]
        
        if length(par_entry) == 1
            p=plot(
            samples, par,
            mean = false, std = false, globalmode = true, marginalmode = true,
            nbins = 200
            ) # TO DO: add a way to constrain the posterior in [0; max from config] or [0; right-est entry on the x axis for signal]
            savefig(joinpath(output, "plots/$(par)_marg_posterior.pdf"))
            ct += 1
            
        # multivariate parameters    
        else
            for idx in 1:length(par_entry) 
                xlab = string("$(par)[$(idx)]")
                ylab = string("P($(par)[$(idx)])")
                
                p=plot(
                unshaped_samples, ct,
                mean = false, std = false, globalmode = true, marginalmode = true,
                nbins = 200, xlabel = xlab, ylabel = ylab,
                )
                savefig(joinpath(output, "plots/$(par)_$(idx)_marg_posterior.pdf"))
                ct += 1
            end
        end
    end
    
    # all parameters together
    plot(
        samples,
        mean = false, std = false, globalmode = false, marginalmode = true,
        nbins = 200
    )
    savefig(joinpath(output, "plots/all_marg_posterior_2D.pdf"))
        
    # fit over data (TO DO)
    """
    # create histo with energies 
    energies = []
    for (idx_k, part_k) in enumerate(partitions[1])
        if events[1][idx_k] != Any[]
            for energy in events[1][idx_k]
                append!(energies, events[1][idx_k])
            end
        end
    end
    hist_data = append!(Histogram(1930:2:2190), energies)
    p_fit = plot_data(hist,name,nothing,samples,fit_function)
    savefig(joinpath(output, "plots/fit_over_data.pdf"))
    """
        
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
    #display(p)
  
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
    #display(p)
  
end