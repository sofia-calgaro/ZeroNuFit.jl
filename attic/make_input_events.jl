using Pkg
Pkg.activate(".")
using JSON
using DataStructures
using PropDicts
using LegendHDF5IO
using TypedTables
using LegendDataManagement

function make_events_file(skmphy_path, file_path, meta_path)
    files =readdir(skmphy_path)
   
    files=[joinpath(skmphy_path,f) for f in files]
    
    vsel=ValiditySelection("20230312T043356Z", :phy)
    chmap = LegendDataManagement.AnyProps(meta_path).hardware.configuration.channelmaps(vsel)
    map = Dict( chmap[name].daq.rawid => String(name) for name in keys(chmap))
    E_fit=[]
    t_fit=[]
    d_fit=[]
    for f in files
        data = lh5open(f)["skm"][:]
        
        data   =data[.!data.coincident.muon_offline .&&
            .!data.coincident.spms .&&
            data.geds.psd.is_bb_like
        ]
    
        E = data.geds.energy
        T=  data.trigger.timestamp
        R= data.geds.rawid
        ranges =[[1930.0, 2098.511], [2108.511, 2113.513], [2123.513, 2190.0]]
       
        for (e,t,ra) in zip(E,T,R)
            for r in ranges
                if (e<r[2] && e>r[1])
                    append!(E_fit,e)
                    append!(t_fit,t)
                    append!(d_fit,[map[Int(ra)]])
                end
               
             end
        end
  
    end

    output=Dict("events" => [])
    for (t,d,e) in zip(t_fit,d_fit,E_fit)
        append!(output["events"],[Dict("experiment"=>"L200","timestamp"=>t,"energy"=>e,"detector"=>d)])
    end
    open(file_path, "w") do file
        write(file, json(output,4))
    end
end

make_events_file("/some_path/generated/tier/skm/phy", "events_example.json","legend-metadata/")
