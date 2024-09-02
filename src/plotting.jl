using Random, LinearAlgebra, Statistics, Distributions, StatsBase
using Plots
using BAT, DensityInterface, IntervalSets
using TypedTables
using Cuba


default(
    titlefont = "Roboto",           # Font for plot titles
    guidefont = "Roboto",               # Font for axis labels
    tickfont = "Arial",              # Font for axis ticks
    legendfont = "Arial",           # Font for the legend
    framestyle=:box,               # Grid line transparency
    background_color = :white   ,       # Background color of the plot,
    titlefontsize=16,     # Global title font size
    guidefontsize=16,     # Global axis label font size
    tickfontsize=10,      # Global tick label font size
    legendfontsize=10     # Global legend font size
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
        title = "$name data",
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
function make_plots(samples,pars,hist,name,fit_function) # ,pdfname=nothing
    samples_mode = mode(samples)
        
    for par in pars
    
        p=plot(
        samples, par,
        mean = false, std = false, globalmode = false, marginalmode = true,
        nbins = 200
        )
        #display(p)
    end
        
    p_fit = plot_data(hist,name,nothing,samples,fit_function,savefig)
    #display(p_fit)
    
    ### TO DO : implement a way to save plots (p, p_fit) in a single PDF!
        
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