"""
Script to print and save some marginalized posterior distributions and compare different background models (if present)
Main Authors: Sofia Calgaro, Toby Dixon
"""
import os,json,h5py,math
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
bin_width_qbb = 0.0005
nbins_qbb = np.arange(xmin, xmax + bin_width_qbb, bin_width_qbb)
bin_width_S = 0.5
xmin_S = np.min(0)
xmax_S = np.max(110)
nbins_S = np.arange(xmin_S, xmax_S + bin_width_S, bin_width_S)
bin_width_S_low = 0.25
nbins_S_low = np.arange(xmin_S, xmax_S + bin_width_S_low, bin_width_S_low)
gm_type = "marginalized_modes"

center = 1930
expo = 42.25
fit_range = [[1930, 2099], [2109, 2114], [2124, 2190]]
Qbb = 2039.06

def normalize():
    range_l = [arr[0] for arr in fit_range] 
    range_h = [arr[1] for arr in fit_range] 
    norm = sum([h - l for h, l in zip(range_h, range_l)])

    return norm

def compute_deltaE():
    deltaE = sum([arr[2] - arr[1] for arr in fit_range])
    return deltaE

def compute_window():
    window = fit_range[-1][-1] - fit_range[0][0]
    return window

def normalize_sq():
    range_l = [arr[0] for arr in fit_range] 
    range_h = [arr[1] for arr in fit_range] 
    window_sq = sum([h**2 - l**2 for h, l in zip(range_h, range_l)])
    return window_sq

def exp_stable(x):
    if abs(x) < 1e-6:
        return 1 + x + (x**2) / 2 + (x**3) / 6
    else:
        return math.exp(x)

norm = normalize()
window = compute_window()
window_sq = normalize_sq()
sum_range = normalize()
sum_range_sq = normalize_sq()
print("norm:", norm)
print("window:", window)

def bkg_at_qbb_linear(bkg_list, slope_list):

    bkg_at_qbb_list = []
    for idx,b in enumerate(bkg_list):
        s = slope_list[idx]
        
        term_1 = (b * expo * norm)
        
        norm_all = sum_range * (1 - s * center / window) + s * sum_range_sq / (2 * window)
        term_2 = (1 + s*(Qbb-center)/window) / norm_all
        
        bkg_at_qbb_list.append(term_1 * term_2)
        
    return bkg_at_qbb_list 



def bkg_at_qbb_uniform(bkg_list):

    bkg_at_qbb_list = []
    for idx,b in enumerate(bkg_list):

        term_1 = (b * expo * norm)
        term_2 =  1 / norm
        
        bkg_at_qbb_list.append(term_1 * term_2)
        
    return bkg_at_qbb_list 


def bkg_at_qbb_exponential(bkg_list, slope_list):
    range_l = [arr[0] for arr in fit_range] 
    range_h = [arr[1] for arr in fit_range]  
    centers = [center, center, center]
    
    bkg_at_qbb_list = []
    for idx,b in enumerate(bkg_list):
        R = slope_list[idx]
        Rt = R / window

        term_1 = (b * expo * norm)
        
        if abs(Rt) > 1e-6:
            norm_all = (-sum([exp_stable(-Rt * (center - l)) for l in range_l]) +
                    sum([exp_stable(-Rt * (center - h)) for h in range_h])) / Rt
        else:
            norm_all = normalize()
        term_2 = exp_stable((Qbb-center)*Rt) / norm_all
        
        bkg_at_qbb_list.append(term_1 * term_2)
    
    return bkg_at_qbb_list 


def smallest_68_ci(bkg, nbins):
    counts, bin_edges = np.histogram(bkg, bins=nbins, density=True)
    max_bin_index = np.argmax(counts)
    max_bin_index = np.argmax(counts)
    max_bin_count = counts[max_bin_index]
    max_bin_start = bin_edges[max_bin_index]
    max_bin_end = bin_edges[max_bin_index + 1]
    marg_mode = (max_bin_end+max_bin_start)/2
    total_area = np.sum(counts * np.diff(bin_edges)) 
    
    target_area = 0.68 * total_area
    
    accumulated_area = 0.0
    left_index = max_bin_index
    right_index = max_bin_index
    bin_widths = np.diff(bin_edges)
    
    while accumulated_area < target_area:
        left_area = counts[left_index] * bin_widths[left_index] if left_index > 0 else 0
        right_area = counts[right_index] * bin_widths[right_index] if right_index < len(counts) - 1 else 0

        if left_area >= right_area and left_index > 0:
            accumulated_area += left_area
            left_index -= 1
        elif right_index < len(counts) - 1:
            accumulated_area += right_area
            right_index += 1
        else:
            break
    
    ci_low = bin_edges[left_index]
    ci_high = bin_edges[right_index + 1]
    
    return marg_mode, ci_low, ci_high
    

def different_bkg_models():
    with PdfPages(f"comparison_of_different_bkg_shapes.pdf") as pdf:
        """
        # l200
        flat_file = f"../ovbb_old/output_old/fit_9_l200_uniform_1BI_CorrEff/mcmc_files/samples.h5"
        linear_file = f"../output/fit_2_l200_linear_CorrEff/mcmc_files/samples.h5"
        expo_file = f"../output/fit_1_l200_expo_CorrEff/mcmc_files/samples.h5"
        bkg_name = "B_l200a_all"
        """
        # gerda
        flat_file = f"../ovbb_old/output_old/wrong_expo_l200/fit_1_gerda_phII_uniform_1BI_CorrEff/mcmc_files/samples.h5"
        linear_file = f"../output/fit_3_gerda_linear_CorrEff/mcmc_files/samples.h5"
        expo_file = f"../output/fit_4_gerda_expo_CorrEff/mcmc_files/samples.h5"
        bkg_name = "B_gerda_all_pII"

        fig, ax = plt.subplots(figsize=(5,3.3))
        signal_flat = read_samples(flat_file, 'S')
        plt.hist(signal_flat, bins=nbins_S, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(signal_flat, bins=nbins_S, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        
        signal_linear = read_samples(linear_file, 'S')
        plt.hist(signal_linear, bins=nbins_S, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(signal_linear, bins=nbins_S, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        
        signal_expo = read_samples(expo_file, 'S')
        plt.hist(signal_expo, bins=nbins_S, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(signal_expo, bins=nbins_S, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.4)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.4)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.4)
        ax.legend(handles=[flt_patch, lin_patch, exp_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlabel(r'S ($10^{-27}$ yr$^{-1}$)')
        plt.ylabel('Probability density')
        plt.title("Comparison of signal posteriors")
        plt.xlim(0,100)
        #plt.yscale("log")
        pdf.savefig(bbox_inches='tight')
        plt.close()
        print("plotted comparison of posterior pdfs for S")
        
        fig, ax = plt.subplots(figsize=(5,3.3))
        plt.hist(signal_flat, bins=nbins_S_low, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(signal_flat, bins=nbins_S_low, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        plt.hist(signal_linear, bins=nbins_S_low, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(signal_linear, bins=nbins_S_low, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        plt.hist(signal_expo, bins=nbins_S_low, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(signal_expo, bins=nbins_S_low, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.4)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.4)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.4)
        ax.legend(handles=[flt_patch, lin_patch, exp_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlabel(r'S ($10^{-27}$ yr$^{-1}$)')
        plt.ylabel('Probability density')
        plt.title("Comparison of signal posteriors - zoom")
        plt.xlim(0,20)
        #plt.yscale("log")
        pdf.savefig(bbox_inches='tight')
        plt.close()
        print("plotted comparison of posterior pdfs for S (zoom)")

        fig, ax = plt.subplots(figsize=(5,3.3))
        bkg_flat = read_samples(flat_file, f'{bkg_name}')
        plt.hist(bkg_flat, bins=nbins, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(bkg_flat, bins=nbins, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        
        bkg_linear = read_samples(linear_file, f'{bkg_name}')
        plt.hist(bkg_linear, bins=nbins, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(bkg_linear, bins=nbins, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        
        bkg_expo = read_samples(expo_file, f'{bkg_name}')
        plt.hist(bkg_expo, bins=nbins, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(bkg_expo, bins=nbins, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.4)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.4)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.4)
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
        #plt.yscale("log")
        pdf.savefig( bbox_inches='tight')
        plt.close()
        print("plotted comparison of posterior pdfs for BI")

        fig, ax = plt.subplots(figsize=(5,3.3))
        bkg_flat = read_samples(flat_file, f'{bkg_name}')
        bkg_flat = bkg_at_qbb_uniform(bkg_flat)
        plt.hist(bkg_flat, bins=nbins_qbb, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(bkg_flat, bins=nbins_qbb, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        marg_mode, ci_low, ci_high = smallest_68_ci(bkg_flat, nbins_qbb)
        print("Marginalized mode (flat bkg):", marg_mode)
        print(f"68% CI: [{ci_low}, {ci_high}]")
        
        bkg_linear = read_samples(linear_file, f'{bkg_name}')
        slope_linear = read_samples(linear_file, f'{bkg_name}_slope')
        bkg_linear = bkg_at_qbb_linear(bkg_linear, slope_linear)
        plt.hist(bkg_linear, bins=nbins_qbb, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(bkg_linear, bins=nbins_qbb, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        marg_mode, ci_low, ci_high = smallest_68_ci(bkg_linear, nbins_qbb)
        print("Marginalized mode (linear bkg):", marg_mode)
        print(f"68% CI: [{ci_low}, {ci_high}]")
        
        bkg_expo = read_samples(expo_file, f'{bkg_name}')
        slope_expo = read_samples(expo_file, f'{bkg_name}_slope')
        bkg_expo = bkg_at_qbb_exponential(bkg_expo, slope_expo)
        plt.hist(bkg_expo, bins=nbins_qbb, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(bkg_expo, bins=nbins_qbb, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        marg_mode, ci_low, ci_high = smallest_68_ci(bkg_expo, nbins_qbb)
        print("Marginalized mode (expo bkg):", marg_mode)
        print(f"68% CI: [{ci_low}, {ci_high}]")
        
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.4)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.4)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.4)
        ax.legend(handles=[flt_patch, lin_patch, exp_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlabel(r'Background at Qbb (cts/keV)')
        plt.title("Comparison of bkg posteriors")
        plt.ylabel('Probability density')
        #plt.yscale("log")
        pdf.savefig( bbox_inches='tight')
        plt.close()
        print("plotted comparison of posterior pdfs for bkg at Qbb")

    

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
        bi_bege = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", f'{bkg_name}')
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
        bi_bege = read_samples(f"ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", f'{bkg_name}')
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
        
        
        
different_bkg_models()