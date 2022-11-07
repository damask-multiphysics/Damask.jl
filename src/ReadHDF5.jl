module ReadHDF5

import HDF5
export read_HDF5, get, read_HDF5_short





function get(fname)
    f = HDF5.h5open(fname,"r")

    get = Dict()
    r = r"increment_([0-9]+)"
    for k in keys(f)
        if occursin(r,k)
            get[k] = Dict()
            
            #Metadata of increments needed? Doesnt work:
            #d = attrs(f[k])
            #meta = NamedTuple{Tuple(Symbol.(keys(d)))}(values(d))
            #get[k]= attach_metadata(get[k],meta)
            #println(typeof(f[k]))
            
            _help_read(f[k], get[k])
            
        end
    end
    
    close(f)
    return get
end
function _help_read(item,get)
    for k in keys(item)
        if item[k] isa HDF5.Group
            get[k]=Dict()
            _help_read(item[k],get[k])
        else
            
            d = attrs(item[k])   
            println(d)            
            meta = NamedTuple{Tuple(Symbol.(keys(d)))}(values(d))
            get[k] = attach_metadata(read(item[k]),meta)
            get[k] = read(item[k])
        
        end    
    end
end

const prefix_inc = "increment_"


struct HDF5_Obj
    file::HDF5.File
    export_setup::Bool#?
    structured::Bool
    cells::Vector{Int32}
    size::Vector{Float32}
    origin::Vector{Float32}
    increments::Vector{String}
    times::Vector{Float32}
    N_constituents::Integer
    N_materialpoints::Integer
    homogenizations::Vector{String}
    phases::Vector{String}
    fields::Vector{String}
    visible::Dict #increments,phases,homogenizations, fields accessable here, so no extra members needed?
    fname::String
    _protected::Bool
end


function read_HDF5(filename::AbstractString)
    #no hdf5 file?
    if ! HDF5.ishdf5(filename)
        error("No HDF5-file") #raise error?
    end
    
    file=HDF5.h5open(filename,"r") # "r": read-only 

    # right version
    version_minor=HDF5.read_attribute(file,"DADF5_version_minor") 
    version_major=HDF5.read_attribute(file, "DADF5_version_major")
    if version_major !=0 || version_minor <12 || version_minor >14
        error("unsupported DADF5 version ",version_major,".",version_minor) 
    end    
    export_setup=true #?
    if version_major == 0 && version_minor < 14
        export_setup = false #?
    end

    
    #Structured
    structured=haskey(HDF5.attributes(file["geometry"]), "cells")
    if structured
        cells=HDF5.read_attribute(file["geometry"],"cells")
        _size=HDF5.read_attribute(file["geometry"],"size") #dont name variable size because conflict
        origin=HDF5.read_attribute(file["geometry"],"origin")
    else
        add_curl=add_divergence=add_gradient=nothing #? 
    end

    #read increments
    increments=[i for i in keys(file) if occursin(Regex("$prefix_inc([0-9]+)"),i)] #alternativ: increments=collect(i for i in keys(file) if occursin(Regex("$prefix([0-9]+)"),i))
    increments=sort!(increments)# sorting algo? not sorted by increment number now

    if length(increments) == 0
        error("no increments found")
    end
    
    times=[round(HDF5.read_attribute(file[i],"t/s")) for i in increments]
    
    (N_constituents,N_materialpoints) = size(file["cell_to/phase"]) #other way around in python, because Julia is column-major

    homogenization_dset = file["cell_to/homogenization"] #Vector{NamedTuple{(:label, :entry), Tuple{String, Int64}}} #TODO astype("str") ?
    homogenization = [homogenization_dset[i][:label] for i in 1:size(homogenization_dset)[1]] #size= (4096,) maybe same as N_materialpoints?
    homogenizations=unique!(sort!(homogenization))# sorting_algo?  first sort, then unique because https://discourse.julialang.org/t/how-can-write-a-function-to-find-unique-elements-in-array-without-any-allocation/34005/4

    phase_dset=file["cell_to/phase"] #size= (1, 4096) why?
    
    phase=[phase_dset[1,i][:label] for i in 1:size(phase_dset)[2]] #maybe N_materialpoints?
    phases=unique!(sort!(phase)) #sorting algo?
    
    fields=String[]
    for c in phases
        val=keys(file[string("/",increments[1],"/phase/",c)]) 
        for v in val
            push!(fields,v)
        end
    end
    for m in homogenizations
        val=keys(file[string("/",increments[1],"/homogenization/",m)]) 
        for v in val
            push!(fields,v)
        end
    end
    fields=unique!(sort!(fields))

    visible =   Dict{String,Vector{String}}("increments"=>    increments,
                "phases" =>         phases,
                "homogenizations" => homogenizations,
                "fields" =>         fields,
    )
    
    #TODO 
    fname=filename
    
    _protected=true
    
    return HDF5_Obj(file
    ,export_setup
    ,structured
    ,cells
    ,_size
    ,origin
    ,increments
    ,times
    ,N_constituents
    ,N_materialpoints
    ,homogenizations
    ,phases
    ,fields
    ,visible
    ,fname
    ,_protected)
end





#******** Testing deletable
struct HDF5_Obj_short
    file::HDF5.File
end

function read_HDF5_short(filename::AbstractString)
    if HDF5.ishdf5(filename)
       return HDF5_Obj_short(HDF5.h5open(filename,"r")) 
    else
        return nothing
    end
end
#********

end   #end of module
