using MoM_Basics, MoM_Kernels, LinearAlgebra
using Test
include("Matrix and Solving.jl")
include("PostProcessing.jl")

SimulationParams.SHOWIMAGE = true
setPrecision!(Float32)

@testset "MoM_Kernels.jl" begin

    # project path
    proj_path = joinpath(@__DIR__, "..")

    @testset "Triangle, RWG" begin
        for IE in [:EFIE, :MFIE, :CFIE]
            inputParameters(;frequency = 20e8, ieT = IE)
            filename = joinpath(proj_path, "meshfiles/Tri.nas")
            meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
            @test true
            ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData; sbfT = :RWG)
            @test true

            ICoeff = test_opt_solving(geosInfo, bfsInfo; source = PlaneWave(π, 0, 0f0, 1f0))
            test_postprocessing(ICoeff, geosInfo)

            ICoeff2 = test_mat_solving(geosInfo, nbf)
            test_postprocessing(ICoeff2, geosInfo)
        end
    end

    @testset "Terahedron, PWC and SWG" begin

        filename = joinpath(proj_path, "meshfiles/Tetra.nas")
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        @testset "PWC, SWG" for vbfT in [:SWG, :PWC]

            ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = vbfT)
            @test true

            setGeosPermittivity!(geosInfo, 2(1-0.001im))
            @test true

            inputParameters(;frequency = 20e8, ieT = :EFIE)
            ICoeff = test_opt_solving(geosInfo, bfsInfo)
            test_postprocessing(ICoeff, geosInfo)

            ICoeff2 = test_mat_solving(geosInfo, nbf)
            test_postprocessing(ICoeff2, geosInfo)

        end
    end

    
    @testset "Hexadron, PWC and RBF" begin

        filename = joinpath(proj_path, "meshfiles/Hexa.nas")
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        @testset "PWC, RBF" for vbfT in [:PWC, :RBF]

            ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = vbfT)
            @test true

            setGeosPermittivity!(geosInfo, 2(1-0.001im))
            @test true

            inputParameters(;frequency = 20e8, ieT = :EFIE)
            ICoeff = test_opt_solving(geosInfo, bfsInfo)
            test_postprocessing(ICoeff, geosInfo)

            ICoeff2 = test_mat_solving(geosInfo, nbf)
            test_postprocessing(ICoeff2, geosInfo)

        end
    end

    @testset "Tetra + Hexadron, PWC" begin

        filename = joinpath(proj_path, "meshfiles/TetraHexa.nas")

        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :PWC)
        @test true

        setGeosPermittivity!(geosInfo, 2(1-0.001im))
        @test true

        inputParameters(;frequency = 20e8, ieT = :EFIE)
        ICoeff = test_opt_solving(geosInfo, bfsInfo)
        test_postprocessing(ICoeff, geosInfo)
        ICoeff2 = test_mat_solving(geosInfo, nbf)
        test_postprocessing(ICoeff2, geosInfo)

    end

    @testset "Tri + Tetra, RWG + SWG, RWG + PWC" begin

        filename = joinpath(proj_path, "meshfiles/TriTetra.nas")

        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        @testset "RWG + SWG, RWG + PWC" for vbfT in [:SWG, :PWC]

            ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, sbfT = :RWG, vbfT = vbfT)
            @test true

            setGeosPermittivity!(geosInfo, 2(1-0.001im))
            @test true

            inputParameters(;frequency = 20e8, ieT = :EFIE)
            ICoeff = test_opt_solving(geosInfo, bfsInfo)
            test_postprocessing(ICoeff, geosInfo)

            ICoeff2 = test_mat_solving(geosInfo, nbf)
            test_postprocessing(ICoeff2, geosInfo)

        end

    end

    @testset "Tri + Hexa, RWG + PWC, RWG + RBF" begin

        filename = joinpath(proj_path, "meshfiles/TriHexa.nas")

        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        @testset "RWG + SWG, RWG + PWC" for vbfT in [:PWC, :RBF]

            ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, sbfT = :RWG, vbfT = vbfT)
            @test true

            setGeosPermittivity!(geosInfo, 2(1-0.001im))
            @test true

            inputParameters(;frequency = 20e8, ieT = :EFIE)
            ICoeff = test_opt_solving(geosInfo, bfsInfo)
            test_postprocessing(ICoeff, geosInfo)

            ICoeff2 = test_mat_solving(geosInfo, nbf)
            test_postprocessing(ICoeff2, geosInfo)

        end

    end

    rm("results"; force = true, recursive = true)
end
