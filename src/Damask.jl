module Damask

import Base: view, get
using HDF5, Metadata, NaturalSort, WriteVTK
export read_HDF5, get, view, view_more, view_less, place, export_VTK


const prefix_inc = "increment_"

"""
    Result

Contains informations and a view of a DADF5 (DAMASK HDF5) file.
Is returned by [`read_HDF5`](@ref)
"""
struct Result
    filename::String
    version_major::Int8
    version_minor::Int8
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
    _protected::Bool
end

"""Presentation of `Result` showing `obj.filename` and `obj.visible`"""
Base.show(io::IO, obj::Result) = begin
    str = "Damask-HDF5-Object:\nFilename:\t$(obj.filename)\nVisible:\tincrements:\t $(obj.visible["increments"])\n\t\tphases:\t\t $(obj.visible["phases"])\n\t\thomogenizations: $(obj.visible["homogenizations"])\n\t\tfields:\t\t $(obj.visible["fields"])\n"
    print(io, str)
end

"""
    view(obj::Result; <keyword arguments>)
    
Return a copy of `obj` with the selected view.
True or "*" as keyword value selects all.
False or String[] as keyword value deselects all.

# Arguments
- `obj::HDF_Obj` : HDF5_obj to create a view of
- `increments::Union{Int,Vector{Int},String,Vector{String},Bool,Nothing}=nothing` : Number(s) of increments to select. Negative Integers allowed eg. -1 for last increment, -2 for second last. Defaults to all.
- `times::Union{<:AbstractFloat,Vector{<:AbstractFloat},String,Vector{String},Bool,Nothing}=nothing,` : Simulation time(s) of increments to select. Defaults to all. Mutually exclusive to increments.
- `phases::Union{String,Vector{String},Bool,Nothing}=nothing,` : Name(s) of phases to select. Defaults to all.
- `homogenizations::Union{String,Vector{String},Bool,Nothing}=nothing,`: Name(s) of homogenizations to select. Defaults to all.
- `fields::Union{String,Vector{String},Bool,Nothing}=nothing,` : Name(s) of fields to select. Defaults to all.
- `protected::Union{Bool,Nothing}=nothing` : Protection status of existing data.

# Examples
Create a view on the last increment
```jldoctest
julia> using Damask
julia> res = read_HDF5("my_file.hdf5")
julia> v = view(res,increments = -1)
```
"""
function view(
    obj::Result;
    increments::Union{Int,Vector{Int},String,Vector{String},Bool,Nothing}=nothing,
    times::Union{<:AbstractFloat,Vector{<:AbstractFloat},String,Vector{String},Bool,Nothing}=nothing,
    phases::Union{String,Vector{String},Bool,Nothing}=nothing,
    homogenizations::Union{String,Vector{String},Bool,Nothing}=nothing,
    fields::Union{String,Vector{String},Bool,Nothing}=nothing,
    protected::Union{Bool,Nothing}=nothing #does nothing now
)
    dup = _manage_view(obj, "set", increments, times, phases, homogenizations, fields)
    if protected !== nothing
        if protected == false
            println("Warning: Modification of existing datasets allowed")
        end
        dup._protected = protected
    end
    return dup
end

"""
    view_more(obj::Result; <keyword arguments>)

Return a copy of `obj` with the selected increased view.
True or "*" as value as keyword value selects all.
False or String[] as keyword value deselects all.

# Arguments
- `obj::HDF_Obj` : HDF5_obj to create a view of
- `increments::Union{Int,Vector{Int},String,Vector{String},Bool,Nothing}=nothing` : Number(s) of increments to select. Negative Integers allowed eg. -1 for last increment. Defaults to all.
- `times::Union{<:AbstractFloat,Vector{<:AbstractFloat},String,Vector{String},Bool,Nothing}=nothing,` : Simulation time(s) of increments to select. Defaults to all. Mutually exclusive to increments.
- `phases::Union{String,Vector{String},Bool,Nothing}=nothing,` : Name(s) of phases to select. Defaults to all.
- `homogenizations::Union{String,Vector{String},Bool,Nothing}=nothing,`: Name(s) of homogenizations to select. Defaults to all.
- `fields::Union{String,Vector{String},Bool,Nothing}=nothing,` : Name(s) of fields to select. Defaults to all.

# Examples
Increasing a view from first increment to increment 0,2,4 :
```jldoctest
julia> using Damask
julia> r = view(read_HDF5("my_file.hdf5"), increments=0)
julia> r2 = view_more(r, increments=[2,4])
```
"""
view_more(
    obj::Result;
    increments::Union{Int,Vector{Int},String,Vector{String},Bool,Nothing}=nothing,
    times::Union{<:AbstractFloat,Vector{<:AbstractFloat},String,Vector{String},Bool,Nothing}=nothing,
    phases::Union{String,Vector{String},Bool,Nothing}=nothing,
    homogenizations::Union{String,Vector{String},Bool,Nothing}=nothing,
    fields::Union{String,Vector{String},Bool,Nothing}=nothing
) = _manage_view(obj, "add", increments, times, phases, homogenizations, fields)

"""
    view_less(obj::Result; <keyword arguments>)
    
Return a copy of `obj` with the selected reduced view.
True or "*" as value as keyword value deselect all.
False or String[] for keyword value does nothing.

# Arguments
- `obj::HDF_Obj` : HDF5_obj to create a view of
- `increments::Union{Int,Vector{Int},String,Vector{String},Bool,Nothing}=nothing` : Number(s) of increments to select. Negative Integers allowed eg. -1 for last increment. Defaults to all.
- `times::Union{<:AbstractFloat,Vector{<:AbstractFloat},String,Vector{String},Bool,Nothing}=nothing,` : Simulation time(s) of increments to select. Defaults to all. Mutually exclusive to increments.
- `phases::Union{String,Vector{String},Bool,Nothing}=nothing,` : Name(s) of phases to select. Defaults to all.
- `homogenizations::Union{String,Vector{String},Bool,Nothing}=nothing,`: Name(s) of homogenizations to select. Defaults to all.
- `fields::Union{String,Vector{String},Bool,Nothing}=nothing,` : Name(s) of fields to select. Defaults to all.

# Examples
Increasing a view from first increment to increment 0,2,4 :
```jldoctest
julia> using Damask
julia> r = view(read_HDF5("my_file.hdf5"), increments=0)
julia> r2 = view_more(r, increments=[2,4])
```
"""
view_less(
    obj::Result;
    increments::Union{Int,Vector{Int},String,Vector{String},Bool,Nothing}=nothing,
    times::Union{<:AbstractFloat,Vector{<:AbstractFloat},String,Vector{String},Bool,Nothing}=nothing,
    phases::Union{String,Vector{String},Bool,Nothing}=nothing,
    homogenizations::Union{String,Vector{String},Bool,Nothing}=nothing,
    fields::Union{String,Vector{String},Bool,Nothing}=nothing
) = _manage_view(obj, "del", increments, times, phases, homogenizations, fields)


"""
    _manage_view(obj::Result; <keyword arguments>)

Create a copy of `obj` and set its view by calling [`_manage_choice`](@ref) on each keyword argument
"""
function _manage_view(
    obj::Result,
    action::String,
    increments::Union{Int,Vector{Int},String,Vector{String},Bool,Nothing}=nothing,
    times::Union{<:AbstractFloat,Vector{<:AbstractFloat},String,Vector{String},Bool,Nothing}=nothing,
    phases::Union{String,Vector{String},Bool,Nothing}=nothing,
    homogenizations::Union{String,Vector{String},Bool,Nothing}=nothing,
    fields::Union{String,Vector{String},Bool,Nothing}=nothing,
)
    if increments !== nothing && times !== nothing
        error("\"increments\" and \"times\" are mutually exclusive")
    end
    dup = deepcopy(obj)
    for (type, input) in [("increments", increments), ("times", times), ("phases", phases), ("homogenizations", homogenizations), ("fields", fields)]
        if input !== nothing
            _manage_choice(dup, input, type, action)
        end
    end
    return dup
end

"""
    _manage_choice(obj::Result, input::Bool, type::String, action::String)

Processes `input` for `type` where type is a keyword argument of [`view`](@ref). Should call [`_set_viewchoice`](@ref) or another [`_manage_choice`]-Method.
"""
function _manage_choice(obj::Result, input::Bool, type::String, action::String)
    if type == "times"
        type = "increments"
    end
    if input == true
        _set_viewchoice(obj, getproperty(obj, Symbol(type)), type, action) #obj.type must have same name
    else
        _set_viewchoice(obj, String[], type, action)
    end
end

function _manage_choice(obj::Result, input::String, type::String, action::String)
    if input == "*"
        _manage_choice(obj, true, type, action)
    elseif type == "times"
        _manage_choice(obj, [parse(Float64, input)], type, action)
    else
        _set_viewchoice(obj, [input], type, action)
    end
end

function _manage_choice(obj::Result, input::Vector{String}, type::String, mode::String)
    if type == "times"
        _manage_choice(obj, [parse(Float64, i) for i in input], type, mode)
    else
        _set_viewchoice(obj, input, type, mode)
    end
end
#only type "increments" possible
function _manage_choice(obj::Result, input::Integer, type::String, action::String)
    _manage_choice(obj, [input], type, action)
end
#only type "increments" possible
function _manage_choice(obj::Result, input::Vector{Int64}, type::String, mode::String)
    choice = String[]
    for i in input
        if i >= 0
            push!(choice, prefix_inc * string(i))
        elseif i < 0 && abs(i) <= length(obj.increments)
            push!(choice, obj.increments[end+i+1])
        end
    end
    _set_viewchoice(obj, choice, type, mode)
end
#only type "times" possible
function _manage_choice(obj::Result, input::Float64, type::String, action::String)
    _manage_choice(obj, [input], type, action)
end
#only type "times" possible
function _manage_choice(obj::Result, input::Vector{<:AbstractFloat}, type::String, action::String)
    choice = String[]
    for (i, time) in enumerate(getproperty(obj, Symbol(type)))
        for j in input
            if isapprox(time, j) #set parameters of isapprox?
                push!(choice, getproperty(obj, Symbol("increments"))[i])
            end
        end
    end
    _set_viewchoice(obj, choice, "increments", action)
end

"""
    _set_viewchoice(obj::Result, choice::Vector{String}, type::String, action::String)

Set view `choice`s of `type` to `obj.visible[type]`. 
`action` must be in ["set", "del", "add"] and specifies if `choice` should be set, added to or deleted from the current view of `obj`. 
"""
function _set_viewchoice(obj::Result, choice::Vector{String}, type::String, action::String)
    existing = Set(obj.visible[type])
    valid = intersect(Set(choice), Set(getproperty(obj, Symbol(type))))
    if action == "set"
        obj.visible[type] = sort!([s for s in valid], lt=natural)#maybe bad: from vector to set to vector
    elseif action == "add"
        add = union(existing, valid)
        obj.visible[type] = sort!([s for s in add], lt=natural)
    elseif action == "del"
        diff = setdiff(existing, valid)
        obj.visible[type] = sort!([s for s in diff], lt=natural)
    end
end

"""
    get(obj::Result, output::Union{String,Vector{String}}="*"; flatten::Bool=true, prune::Bool=true)

Collect data reflecting the group/folder structure in the DADF5 file.
Only datasets in the current view and optionally specified by `output` are collected.

# Arguments
- obj: Result to get data from
- output : Names of the datasets to read.
            Defaults to '*', in which case all datasets are read.
- flatten : Remove singular levels of the folder hierarchy.
            This might be beneficial in case of single increment,
            phase/homogenization, or field. Defaults to True.
- prune : Remove branches with no data. Defaults to True.
"""
function get(
    obj::Result,
    output::Union{String,Vector{String}}="*";
    flatten::Bool=true,
    prune::Bool=true
)
    file = HDF5.h5open(obj.filename, "r")
    dict = Dict()

    all = output == "*" ? true : false
    output = output isa String ? [output] : output

    for inc in obj.visible["increments"]
        dict[inc] = Dict([("phase", Dict()), ("homogenization", Dict()), ("geometry", Dict())])
        for out in keys(file[inc*"/geometry"])
            if all || out in output
                dict[inc]["geometry"][out] = _read(file[inc*"/geometry/"*out])
            end
        end
        for ty in ["phase", "homogenization"]
            for label in obj.visible[ty*"s"]
                dict[inc][ty][label] = Dict()
                for field in keys(file[inc*"/"*ty*"/"*label])
                    if field in obj.visible["fields"]
                        dict[inc][ty][label][field] = Dict()
                        for out in keys(file[inc*"/"*ty*"/"*label*"/"*field])
                            if all || out in output
                                dict[inc][ty][label][field][out] = _read(file[inc*"/"*ty*"/"*label*"/"*field*"/"*out])
                            end
                        end
                    end
                end
            end
        end
    end
    close(file)
    if prune
        dict = _dict_prune(dict)
    end
    if flatten && dict isa Dict
        dict = _dict_flatten(dict)
    end
    return dict
end
"""
    _read(dataset::HDF5.Dataset)

Read a HDF5.Dataset including metadata.
"""
function _read(dataset::HDF5.Dataset)
    d = attrs(dataset)
    meta = NamedTuple{Tuple(Symbol.(keys(d)))}(values(d))
    arr = read(dataset)
    return attach_metadata(arr, meta)
end


"""
    place(obj::Result, output::Union{String,Vector{String}}="*"; constituents::Union{Vector{Int64},Nothing}=nothing, flatten::Bool=true, prune::Bool=true, fill_float::Float64=NaN, fill_int::Int64=0)

Merge data into spatial order that is compatible with the VTK representation.

The returned data structure reflects the group/folder structure in the DADF5 file.
Multi-phase data is fused into a single output.
`place` is equivalent to `get` if only one phase/homogenization
and one constituent is present.

# Arguments
- obj : Result to be worked on
- output : 
    Names of the datasets to read.
    Defaults to '*', in which case all datasets of the current view are exported.
- constituents : 
    Constituents to consider.
    Defaults to nothing, in which case all constituents are considered.
- fill_float :
    Fill value for non-existent entries of floating point type.
    Defaults to NaN.
- fill_int : 
    Fill value for non-existent entries of integer type.
    Defaults to 0.
"""
function place(
    obj::Result,
    output::Union{String,Vector{String}}="*";
    constituents::Union{Vector{Int64},Int64,Nothing}=nothing,
    flatten::Bool=true,
    prune::Bool=true,
    fill_float::Float64=NaN,
    fill_int::Int64=0
)
    all_out = output == "*" ? true : false
    output = output isa String ? [output] : output

    constituents isa Int64 ? [constituents] : constituents
    constituents_ = constituents === nothing ? range(1, obj.N_constituents) : constituents
    if !all(i->(i in 1:obj.N_constituents), constituents_)
        error("constituent number does not exist")
    end

    suffixes = obj.N_constituents == 1 || length(constituents_) < 2 ? [""] : ["#" * string(c) for c in constituents_]

    (at_cell_ph, in_data_ph, at_cell_ho, in_data_ho) = _mappings(obj)
    dict = Dict()
    file = HDF5.h5open(obj.filename, "r")
    for inc in obj.visible["increments"]
        dict[inc] = Dict([("phase", Dict()), ("homogenization", Dict()), ("geometry", Dict())])
        for out in keys(file[inc*"/geometry"])
            if all_out || out in output
                dict[inc]["geometry"][out] = _read(file[inc*"/geometry/"*out])
            end
        end
        for ty in ["phase", "homogenization"]
            for label in obj.visible[ty*"s"]
                for field in keys(file[inc*"/"*ty*"/"*label])
                    if field in obj.visible["fields"]
                        if !(field in keys(dict[inc][ty]))
                            dict[inc][ty][field] = Dict()
                        end
                        for out in keys(file[inc*"/"*ty*"/"*label*"/"*field])
                            if all_out || out in output
                                data = _read(file[inc*"/"*ty*"/"*label*"/"*field*"/"*out])
                                dims = ntuple(i -> :, ndims(data) - 1)
                                type = eltype(parent(data))
                                fill_val = type <: Integer ? Int32(fill_int) : Float32(fill_float)
                                if ty == "phase"
                                    if !(out * suffixes[1] in keys(dict[inc][ty][field]))
                                        for suffix in suffixes
                                            dict[inc][ty][field][out*suffix] = _empty_like(obj, data, fill_val)
                                        end
                                    end
                                    for (c, suffix) in zip(constituents_, suffixes)
                                        dict[inc][ty][field][out*suffix][dims..., at_cell_ph[c][label]] = data[dims..., in_data_ph[c][label]]
                                    end
                                end
                                if ty == "homogenization"
                                    if !(out in keys(dict[inc][ty][field]))
                                        dict[inc][ty][field][out] = _empty_like(obj, data, fill_val)
                                    end
                                    dict[inc][ty][field][out][dims..., at_cell_ho[label]] = data[dims..., in_data_ho[label]]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    close(file)
    if prune
        dict = _dict_prune(dict)
    end
    if flatten && dict isa Dict
        dict = _dict_flatten(dict)
    end
    return dict
end

"""
    export_VTK(obj::Result, output::Union{String,Vector{String}}="*"; mode::String="cell", constituents::Union{Vector{Int64},Nothing}=nothing, fill_float::Float64=NaN, fill_int::Int64=0)

Export to VTK cell data. Only datasets in the current view and optionally specified by `output` are exported.
Create one VTK file (.vti) per visible increment.

# Arguments
- obj : Result to be exported 
- output : 
    Names of the datasets to export to the VTK file.
    Defaults to '*', in which case all datasets of the current view are exported.
- mode : {'cell', 'point'}
    Export in cell format or point format. Only cell format possible now.
    Defaults to 'cell'.
- constituents :  
    optional constituents to consider.
    Defaults to nothing, in which case all constituents are considered.
- fill_float : 
    Fill value for non-existent entries of floating point type.
    Defaults to NaN.
- fill_int :
    Fill value for non-existent entries of integer type.
    Defaults to 0.
"""
function export_VTK(obj::Result,
    output::Union{String,Vector{String}}="*";
    mode::String="cell", #only cell-mode now
    constituents::Union{Vector{Int64},Nothing}=nothing,
    fill_float::Float64=NaN,
    fill_int::Int64=0
)

    if lowercase(mode) == "cell"
        mode = "cell"
    elseif lowercase(mode) == "point"
        mode = "point"
    else
        error("invalid mode: " * mode)
    end

    all = output == "*" ? true : false
    output = output isa String ? [output] : output

    constituents_ = constituents === nothing ? range(1, obj.N_constituents) : constituents
    suffixes = obj.N_constituents == 1 || length(constituents_) < 2 ? [""] : ["#" * string(c) for c in constituents_]

    (at_cell_ph, in_data_ph, at_cell_ho, in_data_ho) = _mappings(obj)

    file = HDF5.h5open(obj.filename, "r")

    cells = read_attribute(file["geometry"], "cells")
    origin = read_attribute(file["geometry"], "origin")
    size = read_attribute(file["geometry"], "size")

    x = (origin[1]:size[1]-origin[1]:cells[1]) / cells[1]
    y = (origin[2]:size[2]-origin[2]:cells[2]) / cells[2]
    z = (origin[3]:size[3]-origin[3]:cells[3]) / cells[3]

    n_digits = ndigits(parse(Int64, split(obj.increments[end], prefix_inc)[end]))

    for inc in obj.visible["increments"]
        k_inc = parse(Int64, split(inc, prefix_inc)[end])
        vtkfile = vtk_grid(split(obj.filename, ".")[1] * "_increm" * string(k_inc, pad=n_digits), x, y, z) #TODO different for mode "point" ?

        vtkfile["created", VTKFieldData()] = read_attribute(file, "creator") * " (" * read_attribute(file, "created") * ")"
        if mode == "cell"
            vtkfile["u"] = _read(file[inc*"/geometry/u_n"])
        else
            vtkfile["u"] = _read(file[inc*"/geometry/u_p"])
        end
        for ty in ["phase", "homogenization"]
            for field in obj.visible["fields"]
                outs = Dict()
                for label in obj.visible[ty*"s"]
                    if field in keys(file[inc*"/"*ty*"/"*label])
                        for out in keys(file[inc*"/"*ty*"/"*label*"/"*field])
                            if all || out in output
                                data = _read(file[inc*"/"*ty*"/"*label*"/"*field*"/"*out])
                                dims = ntuple(i -> :, ndims(data) - 1)
                                fill_val = eltype(parent(data)) <: Integer ? Int32(fill_int) : Float32(fill_float)
                                if ty == "phase"
                                    if !(out * suffixes[1] in keys(outs))
                                        for suffix in suffixes
                                            outs[out*suffix] = _empty_like(obj, data, fill_val)
                                        end
                                    end
                                    for (c, suffix) in zip(constituents_, suffixes)
                                        outs[out*suffix][dims..., at_cell_ph[c][label]] = data[dims..., in_data_ph[c][label]]
                                    end
                                end
                                if ty == "homogenization"
                                    if !(out in keys(outs))
                                        outs[out] = _empty_like(obj, data, fill_val)
                                    end
                                    outs[out][dims..., at_cell_ho[label]] = data[dims..., in_data_ho[label]]
                                end
                            end
                        end
                    end
                end

                for (label, dataset) in outs
                    vtkfile["/"*ty*"/"*field*"/"*label*" / "*dataset.unit] = dataset
                end
            end
        end
        vtk_save(vtkfile)
    end
    close(file)
end

"""
    _empty_like(obj::Result, dataset, fill_val)

Creates `Array` filled with `fill_val` of the dimensions of `dataset` in the first dimensions 
and the amount of `obj.N_materialpoints` in the last dimension. Attach Metadata of the dataset.
"""
function _empty_like(obj::Result, dataset, fill_val)
    shape = (size(dataset)[1:end-1]..., obj.N_materialpoints)
    return attach_metadata(fill(fill_val, shape), metadata(dataset))
end

"""
    _mappings(obj::Result)

Extracts mapping data for phases and homogenizations to place data spatially
"""
function _mappings(obj::Result)
    file = HDF5.h5open(obj.filename, "r")

    at_cell_ph = [Dict{String,Vector{Int64}}() for _ in 1:obj.N_constituents]
    in_data_ph = [Dict{String,Vector{Int64}}() for _ in 1:obj.N_constituents]
    at_cell_ho = Dict{String,Vector{Int64}}()
    in_data_ho = Dict{String,Vector{Int64}}()
    phase_dset = file["cell_to/phase"][1:end, 1:end]
    phase = getindex.(phase_dset, :label)
    homogenization_dset = file["cell_to/homogenization"][1:end]
    homogenization = getindex.(homogenization_dset, :label)

    for c in range(1, obj.N_constituents)
        for label in obj.visible["phases"]
            at_cell_ph[c][label] = findall(x -> x == label, phase[c, :])
            in_data_ph[c][label] = getindex.(phase_dset[c, at_cell_ph[c][label]], :entry) .+ 1 # plus 1 because Julia is 1-based
        end
    end
    for label in obj.visible["homogenizations"]
        at_cell_ho[label] = findall(x -> x == label, homogenization)
        in_data_ho[label] = getindex.(homogenization_dset[at_cell_ho[label]], :entry) .+ 1
    end
    close(file)
    return (at_cell_ph, in_data_ph, at_cell_ho, in_data_ho)
end



"""
    _dict_prune(d::Dict)

Recursively remove empty dictionaries. Return a new pruned `Dict` which may be empty. 
"""
function _dict_prune(d::Dict)
    new = Dict{String,Any}()
    for k in keys(d)
        if d[k] isa Dict
            d[k] = _dict_prune(d[k])
        end
        if !(d[k] isa Dict) || !isempty(d[k])
            new[k] = d[k]
        end
    end
    return new
end

"""
    _dict_flatten(d::Dict)

Recursively remove keys of single-entry dictionaries. Return a new flattened `Dict` or a single Array in `d`.
"""
function _dict_flatten(d::Dict)
    if length(d) == 1
        entry = d[collect(keys(d))[1]]
        new = _dict_flatten(entry)
    else
        new = Dict()
        for k in keys(d)
            new[k] = _dict_flatten(d[k])
        end
    end
    return new
end

_dict_flatten(d::Any) = d

"""
    read_HDF5(filename::String)

Starting point for post-processing. 
Read a DADF5 (DAMASK HDF5) file containing DAMASK results and
get a Result returned needed for further post-processing.

# Examples
```jldoctest
julia> using Damask
julia> res = read_HDF5("my_file.hdf5")
```
If the path/filename contains backslashes use "raw" in front of it or use double backslashes 
```jldoctest
julia> res = read_HDF5(raw"my_file.hdf5")
```
"""
function read_HDF5(filename::String)
    #no hdf5 file?
    if !HDF5.ishdf5(filename)
        throw(Base.IOError("No HDF5-File", 0))
    end
    file = HDF5.h5open(filename, "r") # "r": read-only 

    # right version
    version_minor = Int8(HDF5.read_attribute(file, "DADF5_version_minor"))
    version_major = Int8(HDF5.read_attribute(file, "DADF5_version_major"))
    if version_major != 1 || version_minor != 0
        error("unsupported DADF5 version ", version_major, ".", version_minor)
    end

    #Structured
    structured = haskey(HDF5.attributes(file["geometry"]), "cells")
    if structured
        cells = HDF5.read_attribute(file["geometry"], "cells")
        _size = HDF5.read_attribute(file["geometry"], "size") #dont name variable size because conflict
        origin = HDF5.read_attribute(file["geometry"], "origin")
    end

    #read increments
    increments = [i for i in keys(file) if occursin(Regex("$prefix_inc([0-9]+)"), i)]
    increments = sort!(increments, lt=natural)
    if length(increments) == 0
        error("no increments found")
    end

    times = [round(HDF5.read_attribute(file[i], "t/s")) for i in increments]
    (N_constituents, N_materialpoints) = size(file["cell_to/phase"]) #other way around in python, because Julia is column-major

    homogenization = getindex.(file["cell_to/homogenization"][1:end], :label)
    homogenizations = unique!(sort!(homogenization))

    phase = getindex.(file["cell_to/phase"][1:end, 1:end], :label)
    phases = unique!(sort!(vec(phase))) #sorting algo? vec() for 1D-Array

    fields = String[]
    for c in phases
        append!(fields, keys(file[string("/", increments[1], "/phase/", c)]))
    end
    for m in homogenizations
        append!(fields, keys(file[string("/", increments[1], "/homogenization/", m)]))
    end
    fields = unique!(sort!(fields))
    close(file)
    visible = Dict{String,Vector{String}}("increments" => increments,
        "phases" => phases,
        "homogenizations" => homogenizations,
        "fields" => fields,
    )

    _protected = true

    return Result(filename, version_major, version_minor, structured, cells, _size, origin, increments, times, N_constituents, N_materialpoints, homogenizations, phases, fields, visible, _protected)
end

end   #end of module
