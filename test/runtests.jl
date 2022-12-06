using MoM_Basics, MoM_Kernels
using Test
include("Matrix and Solving.jl")

@testset "MoM_Kernels.jl" begin
    @testset "Triangle, RWG" begin
        filename = "../meshfiles/sphere_r1m_metal_1GHz.nas"
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true
        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData)
        @test true
        for IE in [:EFIE, :MFIE, :CFIE]
            inputParameters(;frequency = 10e8, ieT = IE)
            test_matrix_solving(geosInfo, bfsInfo)
        end
    end

    @testset "Terahedron, PWC and SWG" begin

        filename = "../meshfiles/sphereShellTetra50mm.nas"
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        @testset "PWC, RBF" for vbfT in [:PWC, :SWG]

            ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = vbfT)
            @test true

            setGeosPermittivity!(geosInfo, 2 + 0.001im)
            @test true

            inputParameters(;frequency = 5e8, ieT = :EFIE)
            test_matrix_solving(geosInfo, bfsInfo)

        end
    end

    
    @testset "Hexadron, PWC and RBF" begin

        filename = "../meshfiles/sphereShellHexa50mm.nas"
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        @testset "PWC, RBF" for vbfT in [:PWC, :RBF]

            ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = vbfT)
            @test true

            setGeosPermittivity!(geosInfo, 2 + 0.001im)
            @test true

            inputParameters(;frequency = 5e8, ieT = :EFIE)
            test_matrix_solving(geosInfo, bfsInfo)

        end
    end


    @testset "Tetra + Hexadron, PWC" begin

        filename = "../meshfiles/sphereShellHT50mm.nas"
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :PWC)
        @test true

        setGeosPermittivity!(geosInfo, 2 + 0.001im)
        @test true

        inputParameters(;frequency = 5e8, ieT = :EFIE)
        test_matrix_solving(geosInfo, bfsInfo)

    end

    rm("results"; force = true, recursive = true)
end
