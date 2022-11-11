module ReadHDF5

using HDF5, Metadata, NaturalSort
export read_HDF5, get, view


const prefix_inc = "increment_"


struct HDF5_Obj
    file::String
    export_setup::Bool#?
    structured::Bool
    cells::Vector{Int64}
    size::Vector{Float64}
    origin::Vector{Float32}
    increments::Vector{String}
    times::Vector{Float32}
    N_constituents::Int64
    N_materialpoints::Int64
    homogenizations::Vector{String}
    phases::Vector{String}
    fields::Vector{String}
    visible::Dict{String,Vector{String}} 
    fname::String
    _protected::Bool
end


function view(
    obj::HDF5_Obj;
    action::String="set", #where action in ["set","add","del"], only "set" now
    increments::Union{Int, Vector{Int},String,Vector{String},Bool,Nothing}=nothing,
    times::Union{<:AbstractFloat,Vector{<:AbstractFloat}, String,Vector{String},Bool,Nothing}=nothing, 
    phases::Union{String,Vector{String}, Bool,Nothing}=nothing,
    homogenization::Union{String,Vector{String}, Bool,Nothing}=nothing,
    fields::Union{String,Vector{String}, Bool,Nothing}=nothing,
    protected::Union{Bool,Nothing}=nothing
)
    if increments !== nothing && times !== nothing
        error("\"increments\" and \"times\" are mutually exclusive")
    end

    dup = deepcopy(obj)
    for (type,input) in zip(["increments","times","phases","homogenizations","fields"],[increments,times,phases,homogenization,fields])
        if input !== nothing
           _manage_choice(dup,input,type,"set")  
        end
    end
    if protected !== nothing
        if protected ==false
            println("Warning: Modification of existing datasets allowed")
        end
        dup._protected=protected
    end
    return dup
end

#@enum mode set add del #maybe? or parameter datatype

function _manage_choice(obj::HDF5_Obj, input::Bool,type::String,action::String)
    if type=="times"
        type="increments"
    end
    if input==true
        _set_viewchoice(obj,getproperty(obj,Symbol(type)),type,action)#field name must be the same as type
    else
        _set_viewchoice(obj,String[],type,action) 
    end
end

function _manage_choice(obj::HDF5_Obj, input::String,type::String,action::String)
   if input=="*"
        _manage_choice(obj,true,type,action)
   elseif type=="times"
        _manage_choice(obj,[parse(Float64,i)],type,action)
   else
        _set_viewchoice(obj,[input],type,action)
   end
end
function _manage_choice(obj, input::Vector{String},type::String,mode::String)
    if type=="times"
        _manage_choice(obj,[parse(Float64,i) for i in input],type,mode)
    else
        _set_viewchoice(obj,input,type,mode)
    end
end
#only type "increments" possible
function _manage_choice(obj::HDF5_Obj, input::Int64,type::String,action::String)
    _set_viewchoice(obj,[prefix_inc*string(input)],type,action)
end
#only type "increments" possible
function _manage_choice(obj, input::Vector{Int64},type::String,mode::String)
    _set_viewchoice(obj,[prefix_inc*string(i) for i in input],type,mode)
end
#only type "times" possible
function _manage_choice(obj::HDF5_Obj, input::Float64,type::String,action::String)
    _manage_choice(obj,[input],type,action)
end
#only type "times" possible
function _manage_choice(obj, input::Vector{<:AbstractFloat},type::String,action::String)   
    choice=String[]
    println(getproperty(obj,Symbol(type)))
    for (i,time) in enumerate(getproperty(obj,Symbol(type)))
        for j in input
            if isapprox(time,j) #set parameters of isapprox?
                push!(choice,getproperty(obj,Symbol("increments"))[i])
            end
        end
    end
    _set_viewchoice(obj,choice,"increments",action)
end
function _set_viewchoice(obj::HDF5_Obj, choice::Vector{String}, type::String, action::String)
    existing = Set(obj.visible[type])
    valid=intersect(Set(choice), Set(getproperty(obj,Symbol(type))))
    if action =="set"
        obj.visible[type]=sort!([s for s in valid],lt=natural)      #maybe bad: from vector to set to vector
    elseif action=="add"
        add=union(existing,valid)
        obj.visible[type]=sort!([s for s in valid],lt=natural)
    elseif action =="del"
        diff=setdiff(existing,valid)
        obj.visible[type]=sort!([s for s in valid],lt=natural)
    end
end


function get(obj::HDF5_Obj, output::Union{String, Vector{String}}="*")
    file = HDF5.h5open(obj.fname,"r")
    dict=Dict()
    
    all = output=="*" ? true : false 
    output = output isa String ? [output] : output
    
    for inc in obj.visible["increments"]
        dict[inc]=Dict([("phase",Dict()),("homogenization",Dict()),("geometry",Dict())])
        for out in keys(file[inc*"/geometry"])
            if all || out in output
                r[inc]["geometry"][out]= _read(file[inc*"/geometry/"*out])
            end
        end
        for ty in ["phase","homogenization"]
            for label in obj.visible[ty*"s"]
                dict[inc][ty][label]=Dict()
                for field in keys(file[inc*"/"*ty*"/"*label])
                    if field in obj.visible["fields"]
                        dict[inc][ty][label][field]=Dict()
                        for out in keys(file[inc*"/"*ty*"/"*label*"/"*field])
                            if all || out in output
                                dict[inc][ty][label][field][out]=_read(file[inc*"/"*ty*"/"*label*"/"*field*"/"*out])
                            end
                        end
                    end   
                end
            end
        end
    end
    return dict
end

function _read(dataset)
    d = attrs(dataset)              
    meta = NamedTuple{Tuple(Symbol.(keys(d)))}(values(d))
    return attach_metadata(read(dataset),meta)
end


function read_HDF5(filename::String)
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
    increments=sort!(increments, lt=natural)# sorting algo? not sorted by increment number now

    if length(increments) == 0
        error("no increments found")
    end
    
    times=[round(HDF5.read_attribute(file[i],"t/s")) for i in increments]
    
    
    (N_constituents,N_materialpoints) = size(file["cell_to/phase"]) #other way around in python, because Julia is column-major

    homogenization_dset = file["cell_to/homogenization"] #TODO astype("str") ?
    homogenization = [homogenization_dset[i][:label] for i in 1:size(homogenization_dset)[1]] #size= (4096,) maybe same as N_materialpoints? #always 1D?
    homogenizations=unique!(sort!(homogenization))# sorting_algo?  first sort, then unique because https://discourse.julialang.org/t/how-can-write-a-function-to-find-unique-elements-in-array-without-any-allocation/34005/4

    phase_dset=file["cell_to/phase"] #size= (1, 4096) because can be 2D if more than one phase?
    
    phase=[phase_dset[i,j][:label] for i=1:N_constituents, j=1:N_materialpoints]
    println(size(phase))
    
    phases=unique!(sort!(vec(phase))) #sorting algo? vec() for 1D-Array
    
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
    
    
    fname=filename #TODO from Python
    
    _protected=true
    close(file)
    
    return HDF5_Obj(filename
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

end   #end of module
