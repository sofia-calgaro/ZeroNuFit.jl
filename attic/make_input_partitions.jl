using Pkg
Pkg.activate(".")
using JSON
using DataStructures
using PropDicts
using TypedTables
using LegendDataManagement
using Printf

const N_A = 6.022E23
const m_76 = 75.92E-3 
const sig_units = 1e-27

exposure_tot_no_eff = 60.98+23.51+103.73+64.45
exposure_tot = 13.65+57.49+43.95+34.52
expected_limit = 10/(2.3*m_76/(sig_units*N_A*log(2)*exposure_tot))

function printdb(db)
    if db isa PropDict
        for (key,value) in db
            println(key)
            printdb(value)
        end
    else
        println("\t",db)
        println()
    end
end

function legend_name_to_type(name)
   name=String(name)
    if (string(name[1])=="V")
        return "IC"
    elseif (string(name[1])=="B" )
        return "BG"
    elseif (string(name[1])=="P")
        return "PC"
    elseif (string(name[1])=="C")
        return "SC"
    end
end

function gerda_name_to_type(name)
    if (name[1]*name[2]=="GD")
        return "BG"
    elseif (name[1]*name[2]=="RG" || name[1]*name[2]=="AN")
        return "SC"
    elseif (name[1]*name[2]=="IC")
        return "IC"
    end
end


function make_partitions_file(meta_path,output_path,by_type=false,extra=false,all=false)

    meta = LegendDataManagement.AnyProps(meta_path)
    meta_parts = LegendDataManagement.AnyProps("legend/")

    partitions = meta_parts.datasets.ovbb_partitions_pars_Nu24
    printflag=false
    exposure_tot=0
    exposure_tot_no_eff=0
    livetime_tot=0
    mass_tot=0
    output=Dict()

    if all == false
        if (by_type==false)
            fit_groups=Dict("l200a"=>Dict("range"=>[[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]],"model"=>"uniform",
                "bkg_name"=>"B_l200a_all"))
        else
            fit_groups=Dict()
            if (extra == true)
                for dettype in ["BG","IC","PC","SC"]
                    fit_groups["$(dettype)_l200a"]=Dict("range"=>[[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]],"model"=>"uniform",
                    "bkg_name"=>"B_l200a_$(dettype)")
                end
            else
                for dettype in ["BG","IC","PC"]
                    fit_groups["$(dettype)_l200a"]=Dict("range"=>[[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]],"model"=>"uniform",
                    "bkg_name"=>"B_l200a_$(dettype)")
                end
            end
        end
    else
        fit_groups=Dict()
        for dettype in ["Nu24","extra"]
            fit_groups["$(dettype)_l200a"]=Dict("range"=>[[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]],"model"=>"uniform",
            "bkg_name"=>"B_l200a_$(dettype)")
        end
    end

    expo_nu24 = 0
    expo_extra = 0
    println(fit_groups)

    # loop over detectors
    for group in keys(fit_groups)
        output[group]=[]
        for (detector, detdata) in partitions
            if detector == :default
                continue
            end

            if all == true
                if group=="extra_l200a"
                    if (legend_name_to_type(String(detector))!="SC" && (String(detector))[1:3]!="V05")
                        continue
                    end
                    if String(detector) == "V05266B"
                        continue
                    end
                    if String(detector) == "V05261B"
                        continue
                    end 
                else
                    if (legend_name_to_type(String(detector))=="SC")
                        continue
                    end
                    if ((String(detector) != "V05266B" && String(detector) != "V05261B") && String(detector)[1:3]=="V05")
                        println(detector)
                        continue
                    end
                end
            else
                if (by_type && legend_name_to_type(String(detector))!=group[1]*group[2])
                    continue
                end
            end
            
            # apply defaults
            if partitions.default isa PropDicts.MissingProperty 
                detdata_merge = copy(detdata)
            else
                new = partitions.default
                detdata_merge = merge(new,copy(detdata))
            end

            
            # loop over partitions for this detector
            for (partition, pardata) in detdata_merge
                if partition == :default
                    continue
                end

                if detdata_merge.default isa PropDicts.MissingProperty 
                    pardata_merge = copy(pardata)
                else
                    new = detdata_merge.default
                    pardata_merge = merge(new,copy(pardata))
                end
    
                if printflag == true
                    println("for partition $partition ")
                    println("\n")
                end
                part_temp = OrderedDict("experiment"=>"L200","detector"=> detector,"part_name"=>partition,
                                     "start_ts"=>pardata_merge[:span_in_utc_s][1],
                                     "end_ts"=>pardata_merge[:span_in_utc_s][2])
                
                eff = prod([v[:val] for v in values(pardata_merge[:ovbb_acceptance])])
                eff_sigma = sqrt(sum([v[:unc]/v[:val] for v in values(pardata_merge[:ovbb_acceptance])].^2))    

                enrichment=meta.hardware.detectors.germanium.diodes[detector].production.enrichment.val
                enrichment_err =meta.hardware.detectors.germanium.diodes[detector].production.enrichment.unc
                eff*=enrichment
                eff_sigma=eff*sqrt(eff_sigma^2+(enrichment_err/enrichment)^2)                
                
                exposure_tot += eff * pardata_merge[:livetime_in_s] * meta.hardware.detectors.germanium.diodes[detector].production.mass_in_g/(60*60*24*365.25*1000)
                exposure_tot_no_eff += pardata_merge[:livetime_in_s] * meta.hardware.detectors.germanium.diodes[detector].production.mass_in_g/(60*60*24*365.25*1000)
                if group=="Nu24_l200a"
                    expo_nu24 += pardata_merge[:livetime_in_s] * meta.hardware.detectors.germanium.diodes[detector].production.mass_in_g/(60*60*24*365.25*1000)
                else
                    expo_extra += pardata_merge[:livetime_in_s] * meta.hardware.detectors.germanium.diodes[detector].production.mass_in_g/(60*60*24*365.25*1000)
                end
                livetime_tot += pardata_merge[:livetime_in_s]
                mass_tot += meta.hardware.detectors.germanium.diodes[detector].production.mass_in_g
                part_temp["eff_tot"]=eff
                part_temp["eff_tot_sigma"]=eff_sigma
                part_temp["width"]=pardata_merge[:fwhm_in_keV][:val]/2.355
                part_temp["width_sigma"]=pardata_merge[:fwhm_in_keV][:unc]/2.355
                part_temp["fwhm"]=pardata_merge[:fwhm_in_keV][:val]
                part_temp["fwhm_sigma"]=pardata_merge[:fwhm_in_keV][:unc]
                part_temp["exposure"]=pardata_merge[:livetime_in_s]*meta.hardware.detectors.germanium.diodes[detector].production.mass_in_g/(60*60*24*365.25*1000)
                part_temp["bias"]=pardata_merge[:energy_bias_in_keV][:val]
                part_temp["bias_sigma"]=pardata_merge[:energy_bias_in_keV][:unc]
                part_temp["eff_psd"]=pardata_merge[:ovbb_acceptance][:psd][:val]
                part_temp["eff_psd_sigma"]=pardata_merge[:ovbb_acceptance][:psd][:unc]
                part_temp["eff_qc"]=pardata_merge[:ovbb_acceptance][:quality][:val]
                part_temp["eff_qc_sigma"]=pardata_merge[:ovbb_acceptance][:quality][:unc]
                part_temp["eff_lar"]=partitions[:default][:default][:ovbb_acceptance][:lar][:val]
                part_temp["eff_lar_sigma"]=partitions[:default][:default][:ovbb_acceptance][:lar][:unc]
                part_temp["eff_active_volume"]=detdata_merge[:default][:ovbb_acceptance][:active_volume][:val]
                part_temp["eff_active_volume_sigma"]=detdata_merge[:default][:ovbb_acceptance][:active_volume][:unc]
                part_temp["eff_containment"]=detdata_merge[:default][:ovbb_acceptance][:containment][:val]
                part_temp["eff_containment_sigma"]=detdata_merge[:default][:ovbb_acceptance][:containment][:unc]

                append!(output[group],[part_temp])
    
            end
        end
    end
    expected_limit = 10/(2.3*m_76/(sig_units*N_A*log(2)*exposure_tot))
    println("exposure ", exposure_tot_no_eff)
    println("exposure * efficiency ", exposure_tot)
    println("livetime (s) ", livetime_tot)
    println("mass tot (kg) ", mass_tot)
    println("Expected 2.3 count limit = $expected_limit x 10^26 yr")

    println("Nu24 exposure ", expo_nu24)
    println("Extra exposure ", expo_extra)

    output_full=Dict("fit_groups"=>fit_groups,
                "partitions"=>output)
    open(output_path, "w") do file
        write(file, json(output_full,4))
    end

end

make_partitions_file("legend-metadata", "partitions_example.json", false, true, true)