using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using Plots
using BAT, DensityInterface, IntervalSets
using TypedTables
using Cuba
using Base.Filesystem
using PDFmerger: append_pdf!
using ColorSchemes

default(
    framestyle=:box,               # Grid line transparency
    background_color = :white   ,       # Background color of the plot,
    titlefontsize=10,     # Global title font size
    guidefontsize=10,     # Global axis label font size
    tickfontsize=10,      # Global tick label font size
    legendfontsize=8     # Global legend font size
)
# define some constants -TODO add to config
Qbb = 2039.06 # keV
N_A = 6.022E23
m_76 = 75.92E-3 # kg/mol
deltaE = 240 # keV
sig_units =1e-27 # signal is in units of this

tol_colors=ColorSchemes.tol_muted
color_schemes=Dict(:blue =>[tol_colors[1],tol_colors[3],tol_colors[2]],
                   :default=>BAT.default_colors,
                   :red=>[:red4,:red,:salmon],
                   :green=>[:darkgreen,:chartreuse3,:palegreen1],
                   :purple=>[tol_colors[8],tol_colors[9],tol_colors[7]],
                   :muted=>[:olivedrab,:goldenrod,:indianred1]
)

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
    
    b_name=part_k.bkg_name
    model_s_k = log(2) * N_A * part_k.exposure * (part_k.eff_tot + p.α * part_k.eff_tot_sigma) * (p.S*sig_units) / m_76
    model_b_k = deltaE * part_k.exposure * p[b_name]

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

    ymax = 1.5*maximum(hist.weights)
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
        best_fit_pars = BAT.mode(samples)
        plot!(p,1930:0.1:2190,x -> find_a_name(best_fit_pars,x),label="Fit",lw=2,color="red")
    end
    
    # exclude the gamma lines
    shape_x = [2114,2114,2124,2124]
    shape_x2 = [2099,2099,2109,2109]
    shape_y=[0,ymax,ymax,0]
    
    plot!(p, shape_x, shape_y, fillalpha=0.3,fill=true, line=false,fillcolor=:blue, label=false)
    plot!(p, shape_x2, shape_y, fillalpha=0.3,fill=true, line=false,fillcolor=:blue, label=false)
    
    return p
end



##############################################
##############################################
##############################################
function plot_fit_and_data(partitions, events, part_event_index, samples, pars, output, plotflag, toy_idx=nothing)
    
    # create histo with energies 
    energies = []
    for (idx_k, part_k) in enumerate(partitions)
        if events[idx_k] != Any[]
            for energy in events[idx_k]
                append!(energies, energy)
            end
        end
    end
    hist_data = append!(Histogram(1930:1:2190), energies)
    p_fit = plot_data(hist_data,"",partitions,part_event_index,pars,samples,plotflag)
    
    log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
    savefig(joinpath(output, "plots/fit_over_data$log_suffix.pdf"))
    
end

##############################################
##############################################
##############################################

function plot_correlation_matrix(samples,output;par_names=nothing,toy_idx=nothing)
"""
Plots the correlation matrixs
"""
    unshaped_samples, f_flatten = bat_transform(Vector, samples)
    covariance_matrix = cov(unshaped_samples)
    var = std(unshaped_samples)

    corr =  100*sqrt(covariance_matrix ./ (var .* var'))
    heatmap(corr,  xlabel="Parameter Index", ylabel="Parameter Index", color=:diverging_bwr_40_95_c42_n256,clim=(-100,100),ctitle="Correlation Coefficient")
    
    log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
    savefig(joinpath(output,"plots/correlations$log_suffix.pdf"))
end

function plot_two_dim_posteriors(samples,pars,output;par_names=nothing,toy_idx=nothing)
    first_sample = samples.v[1]

    for par_x in pars
        par_entry = first_sample[par_x]
        
        if length(par_entry) != 1
            continue
        end
        for par_y in pars
            if (par_x==par_y)
                continue
            end
            par_entry = first_sample[par_y]
            if length(par_entry) != 1
                continue
            end
            
            x = get_par_posterior(samples,par_x,idx=nothing)
            y = get_par_posterior(samples,par_y,idx=nothing)

            p=histogram2d(x, y, bins=200, cmap=:batlow, xlabel=par_names[par_x], ylabel=par_names[par_y],
            right_margin = 10Plots.mm)
            savefig(p,"temp.pdf")
            
            log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
            append_pdf!(joinpath(output, "plots/2D_posterior$log_suffix.pdf"), "temp.pdf", cleanup=true)
        end
    end


end


##############################################
##############################################
##############################################

function plot_marginal_distr(partitions,samples,pars,output;sqrt_prior=false,
    priors=nothing,par_names=nothing,plot_config=nothing,s_max=nothing,hier=false,toy_idx=nothing)    
"""
Function to plot 1D and 2D marginalized distributions (and priors)
"""
    
    name = split(output, "output/")[end]
    first_sample = samples.v[1]
    unshaped_samples, f_flatten = bat_transform(Vector, samples)
    @debug "Unshaped samples:", bat_report(unshaped_samples)
    
    # remove old file
    log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
    if isfile(joinpath(output, "plots/marg_posterior$log_suffix.pdf"))
        Filesystem.rm(joinpath(output, "plots/marg_posterior$log_suffix.pdf"),force=true)
    end
    

    # get a color scheme
    if plot_config!=nothing && haskey(plot_config,"scheme")
        color_scheme = color_schemes[Symbol(plot_config["scheme"])] 
    else 
        color_scheme=BAT.default_colors
    end
    if plot_config!=nothing && haskey(plot_config,"alpha")
        alpha=plot_config["alpha"]
    else 
        alpha=1
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

            if (par_names !=nothing)
                xname = par_names[par]
            end

            

            p=plot(
            samples, par,
            mean = false, std = false, globalmode = true, marginalmode = true,
            nbins = 200, xlim=(mini,maximum(post)),colors=color_scheme,alpha=alpha,lw=0,linecolor=:black
            ) 
            xaxis!(xname)
            yaxis!("Probability Density.")
            ylims!(0,ylims()[2])
            
            x=range(mini, stop=maximum(post), length=1000)

            # plot prior
            if priors!=nothing
                if (par==:S && sqrt_prior)
                    y= x-> 1 ./(2*sqrt(s_max*x)) 
                else
                    color="black"
                    if (hier==false)
                        y=pdf(priors[par],x)
                        plot!(x,y,label="prior",color="black")

                    elseif (haskey(priors.pdist,par))

                        y=pdf(priors.pdist[par],x)
                        plot!(x,y,label="prior",color="black")

                    else
                        for i in 1:50
                            rando=rand(priors.pdist)
                            rando = samples.v[i]
                            y=pdf(priors.f(rando)[par],x)
                            color="grey"
                            
                            plot!(x,y,color="grey",alpha=0.3,label=nothing)

                        end
                    end
                end
                
            end

            savefig(p,"temp.pdf")
            log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
            append_pdf!(joinpath(output, "plots/marg_posterior$log_suffix.pdf"), "temp.pdf", cleanup=true)
            ct += 1
            
        # multivariate parameters    
        else
            for idx in 1:length(par_entry) 
                post = get_par_posterior(samples,par,idx=idx)

                xlab = string("$(par)[$(idx)]")
                ylab = string("Probability Density.")
                if (par_names !=nothing)
                    xname = par_names[par][idx]
                end

                p=plot(
                unshaped_samples, ct,
                mean = false, std = false, globalmode = true, marginalmode = true,
                nbins = 200, xlabel = xlab, ylabel = ylab, xlim=(minimum(post),maximum(post)),colors=color_scheme,alpha=alpha,
                linecolor=:black
                )
                mini=minimum(post)
                x=range(mini, stop=maximum(post), length=1000)

                if priors!=nothing
                    
                    if (hier==true)
                        y=pdf(priors.pdist[par].v[idx],x)
                    else                       
                        y=pdf(priors[par].v[idx],x)
                    end
                    plot!(x,y,label="prior",color="black")
                end

                xaxis!(xname)
                ylims!(0,ylims()[2])
                savefig(p,"temp.pdf")
                log_suffix = toy_idx == nothing ? "" : "_$(toy_idx)"
                append_pdf!(joinpath(output, "plots/marg_posterior$log_suffix.pdf"), "temp.pdf", cleanup=true)
                ct += 1
            end
        end
    end
        
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