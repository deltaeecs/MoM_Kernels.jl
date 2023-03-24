
# 近场矩阵分块保存
include("MatrixChunk.jl")
include("ZnearChunks.jl")


ZnearT{CT}  = Union{SparseMatrixCSC{CT, Int}, Transpose{CT, SparseMatrixCSC{CT, Int}}} where {CT}

"""
根据八叉树层信息初始化近场矩阵元
"""
function initialZnearCSC(level, nbf::Int)

    FT  =   typeof(level.cubeEdgel)
    # 本盒子
    cubes   =   level.cubes
    
    # CSC矩阵列起始位置指针
    colPtr      =   ones(Int, nbf + 1)

    # 确定近场元的总数、计算列起始位置指针
    nZnear   =   0
    # 对盒子循环
    for iCube in eachindex(cubes)
        # 盒子信息
        cube    =   cubes[iCube]
        # 盒子内包含的源基函数数目
        nSBF    =   length(cube.bfInterval)
        # 对其邻盒循环
        nTBFNear=   0
        for jNearCube in cube.neighbors
            # 邻盒子信息
            nearCube    =   cubes[jNearCube]
            # 对邻盒子中测试基函数数量累加
            nTBFNear   +=   length(nearCube.bfInterval)
        end # jNearCube
        
        # 写入列起始位置指针
        colPtr[cube.bfInterval] .=   (nZnear + 1) .+ ((0:(nSBF-1)) .* nTBFNear)

        # 累加
        nZnear  +=  nSBF*nTBFNear
    end # iCube
    colPtr[end] =   nZnear + 1

    # 预分配CSC矩阵行索引、矩阵元向量
    rowIndices  =   zeros(Int, nZnear)
    Znear       =   zeros(Complex{FT}, nZnear)

    # 对盒子循环，给行索引、列起始位置指针赋值
    for iCube in eachindex(cubes)
        # 盒子信息
        cube    =   cubes[iCube]

        # 对其邻盒循环计算该盒子的邻盒子共包括的测试基函数索引
        nTBFNear=   0
        nearTBFIndices  =   Int[]
        for jNearCube in cube.neighbors
            # 邻盒子信息
            nearCube    =   cubes[jNearCube]
            # 对邻盒子中测试基函数数量累加
            nTBFNear   +=   length(nearCube.bfInterval)
            # 合并
            nearTBFIndices  =   [nearTBFIndices; nearCube.bfInterval]
        end # jNearCube

        # 对源基函数循环，写入对应的测试基函数行索引
        for inBF in cube.bfInterval
            rowIndices[colPtr[inBF]:(colPtr[inBF] + nTBFNear - 1)]  .=  nearTBFIndices
        end

    end # iCube

    
    return SparseMatrixCSC{Complex{FT}, Int}(nbf, nbf, colPtr, rowIndices, Znear)
end


"""
根据八叉树层信息初始化 CSR 压缩稀疏行 (用 transpose 实现) 近场矩阵元
"""
function initialZnearCSR(level, nbf::Int)
    FT  =   typeof(level.cubeEdgel)
    # 本盒子
    cubes   =   level.cubes
    nCubes  =   length(cubes)
    
    # CSR矩阵行起始位置指针
    rowPtr  =   ones(Int, nbf + 1)

    # 确定近场元的总数、计算列起始位置指针
    nZnear  =   0
    # 对盒子循环
    pmeter  =   Progress(nCubes; desc = "Counting elements number of Znear...", dt = 1)
    for iCube in eachindex(cubes)
        # 盒子信息
        cube    =   cubes[iCube]
        # 盒子内包含的场基函数数目
        nTBF    =   length(cube.bfInterval)
        # 对其邻盒循环
        nSBFNear=   0
        for jNearCube in cube.neighbors
            # 邻盒子信息
            nearCube    =   cubes[jNearCube]
            # 对邻盒子中源基函数数量累加
            nSBFNear   +=   length(nearCube.bfInterval)
        end # jNearCube
        
        # 写入列起始位置指针
        rowPtr[cube.bfInterval] .=   (nZnear + 1) .+ ((0:(nTBF-1)) .* nSBFNear)

        # 累加
        nZnear  +=  nTBF*nSBFNear
        next!(pmeter)
    end # iCube
    rowPtr[end] =   nZnear + 1

    # 预分配 CSR 矩阵行索引、矩阵元向量
    print("分配近场矩阵内存中...")
    colIndices  =   zeros(Int, nZnear)
    Znear       =   zeros(Complex{FT}, nZnear)
    print("完毕！")

    # 对盒子循环，给行索引、列起始位置指针赋值
    for iCube in eachindex(cubes)
        # 盒子信息
        cube    =   cubes[iCube]

        # 对其邻盒循环计算该盒子的邻盒子共包括的测试基函数索引
        nTBFNear=   0
        nearTBFIndices  =   Int[]
        for jNearCube in cube.neighbors
            # 邻盒子信息
            nearCube    =   cubes[jNearCube]
            # 对邻盒子中测试基函数数量累加
            nTBFNear   +=   length(nearCube.bfInterval)
            # 合并
            nearTBFIndices  =   [nearTBFIndices; nearCube.bfInterval]
        end # jNearCube

        # 对源基函数循环，写入对应的测试基函数行索引
        for inBF in cube.bfInterval
            colIndices[rowPtr[inBF]:(rowPtr[inBF] + nTBFNear - 1)]  .=  nearTBFIndices
        end

    end # iCube

    
    return transpose(SparseMatrixCSC{Complex{FT}, Int}(nbf, nbf, rowPtr, colIndices, Znear))
end


"""
根据八叉树盒子信息初始化 cube 对应的近场矩阵元块儿
"""
function initialZnearChunks(cube, cubes::AbstractVector; CT = Complex{Precision.FT})

    # 确定近场元的总数、计算列起始位置指针
    rowIndices  =   cube.bfInterval

    # 本地化盒子
    cubeslw     =   getDArgs2local(cubes)
    # 对其邻盒循环
    colIndices  =   Int[]
    for jNearCube in cube.neighbors
        # 邻盒子信息
        nearCube    =   cubeslw[jNearCube]
        # 对邻盒子中测试基函数数量累加
        append!(colIndices, nearCube.bfInterval)
    end # jNearCube
    unique!(sort!(colIndices))

    return MatrixChunk{CT}(rowIndices, colIndices)
end

"""
计算某一块的矩阵向量乘积
"""
function ZnearChunkMulIVec!(ZnearChunk, resultChunk, IVec)

    copyto!(resultChunk, ZnearChunk * IVec)

    nothing

end

getNUnknown(v::AbstractVector) = length(v)
getNUnknown(vs::AbstractVector{VT}) where {VT<:AbstractVector} = sum(length, vs)

"""
给出参数计算矩阵近场元并保存在CSC矩阵中
"""
function calZnearCSC(level, geosInfo::Vector, 
    bfsInfo::Vector)
    
    println("计算矩阵近场元CSC格式稀疏矩阵中...")
    # 初始化CSC矩阵
    nbf         =   getNUnknown(bfsInfo)
    ZnearCSC    =   begin
        if nthreads() <= 4
            # 线程较小时，采用 CSC 格式不会有问题
            initialZnearCSC(level, nbf) 
        else
            #= 线程较大时，采用 CSC 格式的矩阵乘法速度会被严重
            限制，又Julia无原生的 CSR 故而采用 CSC 转置 =#
            initialZnearCSR(level, nbf)
        end
    end
    @clock "计算近场矩阵" begin
        calZnearCSC!(level, geosInfo, ZnearCSC)
    end
    
    return ZnearCSC
end

"""
根据积分方程类型选择相应
"""
function calZnearCSC!(level, geosInfo::AbstractVector{VSCellT}, 
    ZnearCSC, bfT::Type{BFT} = VSBFTypes.sbfType) where {BFT<:BasisFunctionType, VSCellT<:SurfaceCellType}
    if SimulationParams.ieT == :EFIE
        # 计算 RWG下 的 EFIE 阻抗矩阵
        calZnearCSCEFIE!(level, geosInfo, ZnearCSC, bfT)
    elseif SimulationParams.ieT == :MFIE
        # 计算 RWG下 的 MFIE 阻抗矩阵
        calZnearCSCMFIE!(level, geosInfo, ZnearCSC, bfT)
    elseif SimulationParams.ieT == :CFIE
        # 计算 RWG下 的 CFIE 阻抗矩阵
        calZnearCSCCFIE!(level, geosInfo, ZnearCSC, bfT)
    end
    return nothing
end

function calZnearCSC!(level, geosInfo::AbstractVector{VSCellT}, 
    ZnearCSC, bfT::Type{BFT} = VSBFTypes.vbfType) where {BFT<:BasisFunctionType, VSCellT<:VolumeCellType}
    # 计算 SWG/PWC/RBF 下的 EFIE 阻抗矩阵
    calZnearCSCEFIE!(level, geosInfo, ZnearCSC, bfT)
    nothing
end

function calZnearCSC!(level, geosInfoV::AbstractVector{VT}, 
    ZnearCSC) where {VT<:AbstractVector}
    if eltype(geosInfoV[1]) <: SurfaceCellType

        # 面元、体元
        tris    =   geosInfoV[1]
        geoVs   =   geosInfoV[2]

        # 面积分 - 面积分 部分
        if SimulationParams.ieT == :EFIE
            # 计算 RWG下 的 EFIE 阻抗矩阵
            calZnearCSCEFIE!(level, tris, ZnearCSC, VSBFTypes.sbfType)
        elseif SimulationParams.ieT == :MFIE
            # 计算 RWG下 的 MFIE 阻抗矩阵
            calZnearCSCMFIE!(level, tris, ZnearCSC, VSBFTypes.sbfType)
        elseif SimulationParams.ieT == :CFIE
            # 计算 RWG下 的 CFIE 阻抗矩阵
            calZnearCSCCFIE!(level, tris, ZnearCSC, VSBFTypes.sbfType)
        end

        # 体积分 - 体积分部分
        # 计算 SWG/PWC/RBF 下的 EFIE 阻抗矩阵
        calZnearCSCEFIE!(level, geoVs, ZnearCSC, VSBFTypes.vbfType)

        # 计算 RWG + SWG/PWC/RBF 下的 EFIE 阻抗矩阵
        calZnearCSCEFIE!(level, tris, geoVs, ZnearCSC, VSBFTypes.vbfType)
    else# 没有面源则是体混合网格

        # # 计算各自阻抗矩阵
        # for geosInfo in geosInfoV
        #     calZnearCSCEFIE!(level, geosInfo, ZnearCSC, VSBFTypes.vbfType)
        # end
        # 计算相互之间的阻抗矩阵
        for i in 1:length(geosInfoV)
            geosA = geosInfoV[i]
            for j in (i+1):length(geosInfoV)
                geosB = geosInfoV[j]
                calZnearCSCEFIE!(level, geosA, geosB, ZnearCSC, VSBFTypes.vbfType)
            end
        end
    end
    nothing
end


# 各类信息
include("ZnearEFIE.jl")
include("ZnearMFIE.jl")
include("ZnearCFIE.jl")
include("ZnearEFIEVSIE.jl")
