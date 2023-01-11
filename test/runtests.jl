using Damask
using Test
using JLD2
using SHA



const ref_path = abspath(joinpath(dirname(pathof(Damask)), "..", "test", "resources"))

const filepath1 = joinpath(ref_path, "12grains6x7x8_tensionY.hdf5")
const filepath2 = joinpath(ref_path, "4grains2x4x3_compressionY.hdf5")


@testset "read_HDF5" begin
    @test_throws Base.IOError read_HDF5("nothing.txt")
    t1 = read_HDF5(filepath1)
    t2 = read_HDF5(filepath2)

    properties = [:N_materialpoints, :N_constituents, :cells, :size, :origin, :increments, :times, :phases, :homogenizations, :fields]
    vals1 = [336, 1, [6, 7, 8], [0.75, 0.875, 1.0], [0.0, 0.0, 0.0],
            ["increment_0", "increment_4", "increment_8", "increment_12", "increment_16", "increment_20", "increment_24", "increment_28", "increment_32", "increment_36", "increment_40"],
            [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0], ["pheno_bcc", "pheno_fcc"], ["SX"], ["mechanical"]]
    vals2 = [24, 8, [2, 4, 3], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0],
            ["increment_0", "increment_2", "increment_4", "increment_6", "increment_8", "increment_10"],
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], ["A", "B", "C"], ["RGC"], ["mechanical"]]
    for (test_obj, vals) in [(t1, vals1), (t2, vals2)]
        for (prop, val) in zip(properties, vals)
            @test getproperty(test_obj, prop) == val
        end
    end
end

@testset "_dict_prune" begin
    ad = Dict("inc1" => Dict("a" => Dict()), "inc2" => Dict("a" => [3, 9], "b" => Dict()))
    @test isequal(Damask._dict_prune(ad), Dict("inc2" => Dict("a" => [3, 9])))
end


@testset "_dict_flatten" begin
    ad = Dict("inc1" => Dict("a" => Dict("b" => [3, 45])), "inc2" => Dict("a" => [3, 9]))
    @test isequal(Damask._dict_flatten(ad), Dict("inc1" => [3, 45], "inc2" => [3, 9]))
end

@testset "view, view_more, view_less" begin
    all_visible_obj = read_HDF5(filepath1)
    t = read_HDF5(filepath1)
    for test_obj in [view(t, times=20.0, fields=false, homogenizations=true), view(t, times=false, phases=true, fields=String[])]

        @testset "view_all" begin

            for type in [:increments, :phases, :fields, :homogenizations]
                #test true and "*" as Parameters for view, view_more, view_less
                for what in [true, "*"]
                    kwargs = Dict()
                    kwargs[type] = what
                    @test isequal(all_visible_obj.visible[String(type)], (view(test_obj; kwargs...)).visible[String(type)])
                    @test isequal(all_visible_obj.visible[String(type)], (view_more(test_obj; kwargs...)).visible[String(type)])
                    @test isequal([], view_less(test_obj; kwargs...).visible[String(type)])
                    #test by setting property List of Objekt as view-Parameter
                    kwargs[type] = getproperty(test_obj, type)
                    @test isequal(all_visible_obj.visible[String(type)], view(test_obj; kwargs...).visible[String(type)])
                    @test isequal(all_visible_obj.visible[String(type)], view_more(test_obj; kwargs...).visible[String(type)])
                    @test isequal([], view_less(test_obj; kwargs...).visible[String(type)])
                end
            end
            #times tested seperatly
            @test isequal(test_obj.increments, view(test_obj, times=test_obj.times).visible["increments"])
            @test isequal(test_obj.increments, view(test_obj, times=true).visible["increments"])
            @test isequal(test_obj.increments, view(test_obj, times="*").visible["increments"])
        end

        @testset "view_none" begin
            #test false and String[] as Parameters
            for type in [:increments, :times, :phases, :fields, :homogenizations], what in [false, String[]]
                kwargs = Dict()
                kwargs[type] = what
                type = type == :times ? :increments : type
                @test isequal([], view(test_obj; kwargs...).visible[String(type)])
            end
        end
    end
    @testset "view_some" begin
        @test isequal(["increment_4"], view(all_visible_obj, times="2.0", phases=false).visible["increments"])
        @test isequal(["increment_4", "increment_28"], view(all_visible_obj, times=["2.0", "14"], phases=String[]).visible["increments"])
        @test isequal(["increment_4"], view(all_visible_obj, times=2.0, phases="pheno_bcc").visible["increments"])
        @test isequal(["increment_4", "increment_28"], view(all_visible_obj, times=[2.0, 14.0], phases=true).visible["increments"])
        @test isequal(["increment_4", "increment_28"], view(all_visible_obj, times=[2.0, 14.0]).visible["increments"])
        @test isequal(["increment_4"], view(all_visible_obj, increments=4).visible["increments"])
        @test isequal(["increment_4", "increment_28"], view(all_visible_obj, increments=[4, 28]).visible["increments"])
        @test isequal(["pheno_fcc"], view_less(all_visible_obj, phases="pheno_bcc", homogenizations=false).visible["phases"])
        @test isequal(["pheno_bcc", "pheno_fcc"], view_more(all_visible_obj, phases="pheno_bcc", homogenizations=false).visible["phases"])
        @test isequal(["increment_40"], view(all_visible_obj, increments=-1).visible["increments"])
        @test isequal(["increment_36", "increment_40"], view(all_visible_obj, increments=[-1, -2]).visible["increments"])
        @test isequal(["increment_40"], view(all_visible_obj, increments=[-1, -30]).visible["increments"])

    end
end

@testset "get" begin
    obj = read_HDF5(filepath2)
    params_view = [Dict(),
        Dict(:increments => 3),
        Dict(:increments => [1, 8, 3, 4, 5, 6, 7]),
        Dict(:phases => ["A", "B"]),
        Dict(:phases => ["A", "C"], :homogenizations => false),
        Dict(:phases => false, :homogenizations => false),
        Dict(:phases => false),
        Dict()]
    #what = [:output, :flatten, :prune]
    params_get = [(["F", "P", "F", "L_p", "F_e", "F_p"], true, true),
        ("F", true, true),
        (["F", "P"], true, true),
        (["F", "P"], true, true),
        (["F", "P", "O"], true, true),
        (["F", "P", "O"], true, true),
        (["Delta_V"], true, true),
        (["u_p", "u_n"], false, false)]
    for i in 1:8
        t_obj = view(obj; params_view[i]...)
        get_obj = get(t_obj, params_get[i][1], flatten=params_get[i][2], prune=params_get[i][3])
        path = joinpath(ref_path, "get", "test_get[$i].jld2")

        #save(path, get_obj, compress=true)
        ref = load(path)
        @test isequal(ref, get_obj)
    end
end

@testset "place" begin
    obj = read_HDF5(filepath2)
    params_view = [Dict(),
        Dict(:increments => 3),
        Dict(:increments => [1, 8, 3, 4, 5, 6, 7]),
        Dict(:phases => ["A", "B"]),
        Dict(:phases => ["A", "C"], :homogenizations => false),
        Dict(:phases => false, :homogenizations => false),
        Dict(:phases => false),
        Dict()]
    #what = [:output, :flatten, :prune, :constituents]
    params_place = [(["F", "P", "F", "L_p", "F_e", "F_p"], true, true, nothing),
        ("F", true, true, [1, 2, 3, 4, 5, 6, 7, 8]),
        (["F", "P"], true, true, 1),
        (["F", "P"], true, true, [1, 2]),
        (["F", "P", "O"], true, true, [1, 8]),
        (["F", "P", "O"], true, true, [1, 2, 3, 4]),
        (["Delta_V"], true, true, [1, 2, 4]),
        (["u_p", "u_n"], false, false, nothing)]
    for i in 1:8
        t_obj = view(obj; params_view[i]...)
        place_obj = place(t_obj, params_place[i][1], flatten=params_place[i][2], prune=params_place[i][3],constituents=params_place[i][4])
        path = joinpath(ref_path, "place", "test_place[$i].jld2")
        #save(path, place_obj, compress=true)
        ref = load(path)
        @test isequal(ref, place_obj)
    end
end

function create_VTK_reference()
    r = view(Damask.read_HDF5(filepath1), increments=40)
    Damask.export_VTK(r)
    val = bytes2hex(sha2_256(read(joinpath(dirname(pathof(Damask)), "..", "12grains6x7x8_tensionY_increm40.vti"))))
    write(joinpath(ref_path, "export_VTK", "12grains6x7x8_tensionY_increm40.vti.sha256"), val)
end

#create_VTK_reference()

@testset "export_VTK" begin
    tmp_path = mktempdir()
    cp(joinpath(ref_path, "12grains6x7x8_tensionY.hdf5"), joinpath(tmp_path, "12grains6x7x8_tensionY.hdf5"))
    r = Damask.read_HDF5(joinpath(tmp_path, "12grains6x7x8_tensionY.hdf5"))
    Damask.export_VTK(r)
    current = bytes2hex(sha2_256(read(joinpath(tmp_path, "12grains6x7x8_tensionY_increm40.vti"), String)))
    reference = read(joinpath(ref_path, "export_VTK", "12grains6x7x8_tensionY_increm40.vti.sha256"), String)
    @test reference == current
end
