"""
Script to print and save some marginalized posterior distributions and compare different background models (if present)
Main Authors: Sofia Calgaro, Toby Dixon
"""
import os,json,h5py
import numpy as np
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

### some colours
c_pistachio = (0.58, 0.87, 0.45)
c_columbiablue = (0.61, 0.87, 1.0)
c_dark_columbiablue = (0.61, 0.87*0.85, 1)
c_frenchblue = (0.0, 0.45, 0.73)
c_dark_gray = (0.33, 0.33, 0.33)
c_dark_pistachio = (0.58*0.85, 0.87*0.85, 0.45*0.85)
c_lava = (0.94, 0.01, 0.05)

### function to read output h5 files with saved posteriors
def read_samples(path,par):
    out=None
    with h5py.File(path, 'r') as hdf:
        
        v = hdf[f"v/{par}"][:]
        w = hdf[f"weight"][:]
        
        out= np.repeat(v, w)
    return out

bin_width = 0.000025
xmin = np.min(0)
xmax = np.max(0.1)
nbins = np.arange(xmin, xmax + bin_width, bin_width)
bin_width_S = 0.5
xmin_S = np.min(0)
xmax_S = np.max(110)
nbins_S = np.arange(xmin_S, xmax_S + bin_width_S, bin_width_S)
gm_type = "marginalized_modes"


def different_bkg_models():
    with PdfPages(f"comparison_of_different_bkg_shapes.pdf") as pdf:
        
        fig, ax = plt.subplots(figsize=(5,3.3))
        signal_flat = read_samples(f"ZeroNuFit.jl/output/fit_9_l200_uniform_1BI_CorrEff/mcmc_files/samples.h5", 'S')
        plt.hist(signal_flat, bins=nbins_S, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(signal_flat, bins=nbins_S, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        signal_linear = read_samples(f"../ZeroNuFit.jl/output/fit_2_l200_linear_CorrEff/mcmc_files/samples.h5", 'S')
        plt.hist(signal_linear, bins=nbins_S, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(signal_linear, bins=nbins_S, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        signal_expo = read_samples(f"../ZeroNuFit.jl/output/fit_1_l200_expo_CorrEff/mcmc_files/samples.h5", 'S')
        plt.hist(signal_expo, bins=nbins_S, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(signal_expo, bins=nbins_S, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.2)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.2)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.2)
        ax.legend(handles=[flt_patch, lin_patch, exp_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlabel(r'S ($10^{-27}$ yr$^{-1}$)')
        plt.ylabel('Probability density')
        plt.title("Comparison of signal posteriors")
        plt.xlim(0,100)
        pdf.savefig(bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(figsize=(5,3.3))
        signal_flat = read_samples(f"ZeroNuFit.jl/output/fit_9_l200_uniform_1BI_CorrEff/mcmc_files/samples.h5", 'B_l200a_all')
        plt.hist(signal_flat, bins=nbins, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(signal_flat, bins=nbins, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        signal_linear = read_samples(f"../ZeroNuFit.jl/output/fit_2_l200_linear_CorrEff/mcmc_files/samples.h5", 'B_l200a_all')
        plt.hist(signal_linear, bins=nbins, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(signal_linear, bins=nbins, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        signal_expo = read_samples(f"../ZeroNuFit.jl/output/fit_1_l200_expo_CorrEff/mcmc_files/samples.h5", 'B_l200a_all')
        plt.hist(signal_expo, bins=nbins, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(signal_expo, bins=nbins, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.2)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.2)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.2)
        ax.legend(handles=[flt_patch, lin_patch, exp_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlim(0,30e-4)
        x_ticks = ax.get_xticks()
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        ax.set_xticks(np.arange(0, 30e-4, 10e-4))
        ax.set_xticklabels(np.arange(0, 30e-4, 10e-4).astype(int))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        plt.xlabel(r'BI ($10^{-4}$ ckky)')
        plt.title("Comparison of BI posteriors")
        plt.ylabel('Probability density')
        pdf.savefig( bbox_inches='tight')
        plt.close()

    

def sig_and_bkg():
    exp = "fit_13_gerdaIandII_5BI_l200_1BI_mjd_uniform_1BI_CorrEffjson_v2" # change me!!
    bkg_name = 'B_mjd-mod2' ### change me!!
    
    with PdfPages(f"signal_and_bkg_posteriors.pdf") as pdf:
        
        fig, ax = plt.subplots(figsize=(5,3.3))
        json_file = json.load(open(f"ZeroNuFit.jl/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file['refined_global_modes']
        c68 = json_file['ci_68']
        bi_all = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", bkg_name)
        plt.hist(bi_all, bins=nbins, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-12)
        plt.hist(bi_all, bins=nbins, histtype='step', density=True, color='darkorange', lw=1, zorder=-12)
        plt.axvline(gms[bkg_name], color='navy', linestyle="-", label="Global mode")
        counts, bin_edges = np.histogram(bi_all, bins=nbins, density=True)
        max_bin_index = np.argmax(counts)
        max_bin_count = counts[max_bin_index]
        max_bin_start = bin_edges[max_bin_index]
        max_bin_end = bin_edges[max_bin_index + 1]
        print("Marginalized mode:", (max_bin_end+max_bin_start)/2*1e4)
        plt.axvline((max_bin_end+max_bin_start)/2, color='r', linestyle="--", label="Marginalized mode")
        plt.axvline(c68[bkg_name][0]['left'], color='r', linestyle=":", label="Smallest 68% CI")
        plt.axvline(c68[bkg_name][0]['right'], color='r', linestyle=":")
        plt.xlim(0,20e-4)
        x_ticks = ax.get_xticks()
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        ax.set_xticks(np.arange(0, 20e-4, 5e-4))
        ax.set_xticklabels(np.arange(0, 20e-4, 5e-4).astype(int))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        plt.xlabel(r'BI ($10^{-4}$ ckky)')
        plt.ylabel('Probability density')
        plt.legend(loc="upper right", frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()


        fig, ax = plt.subplots(figsize=(5,3.3))
        gms = json_file['refined_global_modes']
        c68 = json_file['ci_68']
        s_all = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'S')
        plt.hist(s_all, bins=nbins_S, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-12)
        plt.hist(s_all, bins=nbins_S, histtype='step', density=True, color='darkorange', lw=1, zorder=-12)
        plt.axvline(gms['S'], color='navy', linestyle="-", label="Global mode")
        counts, bin_edges = np.histogram(s_all, bins=nbins_S, density=True)
        max_bin_index = np.argmax(counts)
        max_bin_count = counts[max_bin_index]
        max_bin_start = bin_edges[max_bin_index]
        max_bin_end = bin_edges[max_bin_index + 1]
        plt.axvline((max_bin_end+max_bin_start)/2, color='r', linestyle="--", label="Marginalized mode")
        plt.axvline(c68['S'][0]['left'], color='r', linestyle=":", label="Smallest 68% CI")
        plt.axvline(c68['S'][0]['right'], color='r', linestyle=":")
        plt.xlim(0,70)
        plt.xlabel(r'S ($10^{-27}$ yr$^{-1}$)')
        plt.ylabel('Probability density')
        plt.legend(loc="upper right", frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()


def comparison_bkg():
    
    with PdfPages(f"comparison_of_different_BI_models.pdf") as pdf:
        
        #################### 3 UNCORR  BI ####################
        exp = "fit_10_l200_uniform_3BI_CorrEff" # change me!!
        fig, ax = plt.subplots(figsize=(5,3.3))
        json_file = json.load(open(f"ZeroNuFit.jl/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file[gm_type]
        bi_bege = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_BG')
        plt.hist(bi_bege, bins=nbins, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(bi_bege, bins=nbins, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        bi_icpc = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_IC')
        plt.hist(bi_icpc, bins=nbins, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(bi_icpc, bins=nbins, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        bi_ppc = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_PC')
        plt.hist(bi_ppc, bins=nbins, histtype='stepfilled', density=True, color=c_lava, alpha=0.15, zorder=-6)
        plt.hist(bi_ppc, bins=nbins, histtype='step', density=True, color=c_lava, lw=1, zorder=-6)
        #################### TOTAL ####################
        exp = "fit_9_l200_uniform_1BI_CorrEff" # change me!!
        json_file = json.load(open(f"ZeroNuFit.jl/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file[gm_type]
        bi_bege = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_all')
        plt.hist(bi_bege, bins=nbins, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-12)
        plt.hist(bi_bege, bins=nbins, histtype='step', density=True, color='darkorange', lw=1, zorder=-12)
        icpc_patch = mpatches.Patch(color=c_dark_pistachio, label='Uncorr. BI, ICPC (5 evts)', alpha=0.2)
        #################### style ####################
        bege_patch = mpatches.Patch(color=c_frenchblue, label='Uncorr. BI, BEGe (1 evts)', alpha=0.2)
        ppc_patch = mpatches.Patch(color=c_lava, label='Uncorr. BI, PPC (1 evts)', alpha=0.2)
        all_patch = mpatches.Patch(color='darkorange', label='Single BI (7 evts)', alpha=0.2)
        ax.legend(handles=[all_patch, bege_patch, icpc_patch, ppc_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlim(0,50e-4)
        x_ticks = ax.get_xticks()
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        ax.set_xticks(np.arange(0, 50e-4, 10e-4))
        ax.set_xticklabels(np.arange(0, 50e-4, 10e-4).astype(int))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        plt.xlabel(r'BI ($10^{-4}$ ckky)')
        plt.ylabel('Probability density')
        pdf.savefig(bbox_inches='tight')
        plt.close()


        #################### 3 CORR BI ####################
        exp = "fit_14_l200_uniform_3BIhier_CorrEff" # change me!!
        fig, ax = plt.subplots(figsize=(5,3.3))
        json_file = json.load(open(f"ZeroNuFit.jl/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file[gm_type]
        bi_bege = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_BG')
        plt.hist(bi_bege, bins=nbins, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(bi_bege, bins=nbins, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)

        bi_icpc = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_IC')
        plt.hist(bi_icpc, bins=nbins, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(bi_icpc, bins=nbins, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)

        bi_ppc = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_PC')
        plt.hist(bi_ppc, bins=nbins, histtype='stepfilled', density=True, color=c_lava, alpha=0.15, zorder=-6)
        plt.hist(bi_ppc, bins=nbins, histtype='step', density=True, color=c_lava, lw=1, zorder=-6)
        #################### TOTAL ####################
        exp = "fit_9_l200_uniform_1BI_CorrEff" # change me!!
        json_file = json.load(open(f"ZeroNuFit.jl/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file[gm_type]
        bi_bege = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_all')
        plt.hist(bi_bege, bins=nbins, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-12)
        plt.hist(bi_bege, bins=nbins, histtype='step', density=True, color='darkorange', lw=1, zorder=-12)
        icpc_patch = mpatches.Patch(color=c_dark_pistachio, label='Corr. BI, ICPC (5 evts)', alpha=0.2)
        bege_patch = mpatches.Patch(color=c_frenchblue, label='Corr. BI, BEGe (1 evts)', alpha=0.2)
        ppc_patch = mpatches.Patch(color=c_lava, label='Corr. BI, PPC (1 evts)', alpha=0.2)
        all_patch = mpatches.Patch(color='darkorange', label='Single BI (7 evts)', alpha=0.2)
        ax.legend(handles=[all_patch, bege_patch, icpc_patch, ppc_patch], loc='upper right', ncol=1, frameon=False)
        #################### style ####################
        plt.xlim(0,50e-4)
        x_ticks = ax.get_xticks()
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        ax.set_xticks(np.arange(0, 50e-4, 10e-4))
        ax.set_xticklabels(np.arange(0, 50e-4, 10e-4).astype(int))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        plt.xlabel(r'BI ($10^{-4}$ ckky)')
        plt.ylabel(f'Counts / ({bin_width*1e4}' + r'$\times 10^{-4}$ ckky)')
        plt.ylabel('Probability density')
        pdf.savefig(bbox_inches='tight')
        plt.close()