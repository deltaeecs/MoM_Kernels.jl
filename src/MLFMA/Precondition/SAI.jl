# 本文件提供计算稀疏近似逆预条件的函数

struct SAIPrec{T} <: AbstractMatrix{T}
    mat::AbstractMatrix{T}
end

Base.show(io::IO, unused::MIME{Symbol("text/plain")}, M::SAIPrec{T}) where {T} = show(io, unused, M.mat)
Base.display(M::SAIPrec{T}) where {T} = display(M.mat)
Base.size(M::SAIPrec{T}) where {T} = size(M.mat)

SparseArrays.nnz(M::SAIPrec{T}) where {T} = nnz(M.mat)


"""
x .= M * x
"""
function LinearAlgebra.ldiv!(M::SAIPrec{T}, x::AbstractVector) where {T}
    x .= M.mat * x
end


"""
y .= M * x
"""
function LinearAlgebra.ldiv!(y::AbstractVector, M::SAIPrec{T}, x::AbstractVector) where {T}
    mul!(y, M.mat, x)
end

Base.eltype(::SAIPrec{T}) where {T} = T

Base.:\(M::SAIPrec{T}, x::AbstractVector) where {T}= ldiv!(copy(x), M, x)

"""
计算稀疏近似逆 (Sparse Approximate Inverse (SAI)) 的函数
输入为近场阻抗矩阵CSC, 叶层信息 (也可以为非叶层, 但计算量更大) 
ZnearCSC::ZnearT{CT}
cubes
该函数提供左预条件
"""
function sparseApproximateInversePl(ZnearCSC::ZnearT{CT}, cubes::AbstractVector) where {FT<:Real, CT<:Complex{FT}}
    # 将本函数内的BLAS设为单线程
    nthds = nthreads()
    BLAS.set_num_threads(1)
    # 首先预分配结果, 采用与 ZnearCSC 相同的稀疏模式
    preM    =   deepcopy(ZnearCSC)

    #TODO: 测试使用 CSR 格式矩阵做左预条件会不会快一点
    
    # 所有的盒子
    nCubes  =   length(cubes)

    # 进度条
    pmeter  =   Progress(nCubes, "Pₗ")

    # 对所有盒子循环
    @threads for iCube in 1:nCubes

        # 本盒子与所有邻盒子id
        cube    =   cubes[iCube]
        ineiIDs =   cube.neighbors

        # 本盒子所有基函数的数量
        nNeibfs =   0
        for ineiID in ineiIDs
            nNeibfs    +=  length(cubes[ineiID].bfInterval)
        end
        # 本盒子所有基函数的 id
        neibfs  =   Vector{Int}(undef, nNeibfs)
        nNeibfs   =   0
        for ineiID in ineiIDs
            n   =   length(cubes[ineiID].bfInterval)
            neibfs[(nNeibfs + 1):(nNeibfs + n)]    .=  cubes[ineiID].bfInterval
            nNeibfs   +=  n
        end
        # 本盒子所有基函数的 id
        # cneibfs::Vector{Int}     =   vcat([collect(cubes[ineiID].bfInterval) for ineiID in ineiIDs]...)
        
        # 本盒子与所有邻盒子及其邻盒子的数量
        niNeiNeis=   0
        for ineiID in ineiIDs
            niNeiNeis    +=  length(cubes[ineiID].neighbors)
        end
        # 本盒子与所有邻盒子 id
        iNeisNeiIDs =   Vector{Int}(undef, niNeiNeis)
        niNeiNeis   =   0
        for ineiID in ineiIDs
            n   =   length(cubes[ineiID].neighbors)
            iNeisNeiIDs[(niNeiNeis + 1):(niNeiNeis + n)]    .=  cubes[ineiID].neighbors
            niNeiNeis   +=  n
        end
        unique!(sort!(iNeisNeiIDs))
        # iNeisNeiIDs::Vector{Int} =   unique(sort!(vcat([cubes[ineiID].neighbors for ineiID in ineiIDs]...)))

        # 本盒子所有基函数的数量
        nNeisNeibfs =   0
        for ineiID in iNeisNeiIDs
            nNeisNeibfs    +=  length(cubes[ineiID].bfInterval)
        end
        # 本盒子所有基函数的 id
        neisNeibfs  =   Vector{Int}(undef, nNeisNeibfs)
        nNeisNeibfs   =   0
        for ineiID in iNeisNeiIDs
            n   =   length(cubes[ineiID].bfInterval)
            neisNeibfs[(nNeisNeibfs + 1):(nNeisNeibfs + n)]    .=  cubes[ineiID].bfInterval
            nNeisNeibfs   +=  n
        end
        
        # 本盒子以及邻盒子的基函数 ids 数量
        # cnneibfs::Vector{Int}    =   vcat([collect(cubes[ineiID].bfInterval) for ineiID in iNeisNeiIDs]...)
        
        # 本盒子基函数
        cbfs        =   cube.bfInterval
        # 本盒子基函数起始点在 nneibfs 的位置
        cbfsInCnnei =   cbfs .+ (searchsortedfirst(neisNeibfs, cbfs.start) - cbfs.start)

        # 提取对应的阻抗矩阵
        Znn     =   zeros(CT, nNeibfs, nNeisNeibfs)
        for j in 1:length(neisNeibfs), i in 1:length(neibfs)
            Znn[i, j]  =   ZnearCSC[neibfs[i], neisNeibfs[j]]
        end
        ZnnH    =   Znn'
        
        # cbfsInCnnei =   cbfs .+ (searchsortedfirst(neibfs, cbfs.start) - cbfs.start)

        # Znn     =   zeros(CT, nNeibfs, nNeibfs)
        # for j in 1:length(neibfs), i in 1:length(neibfs)
        #     Znn[i, j]  =   ZnearCSC[neibfs[i], neibfs[j]]
        # end
        # ZnnH    =   Znn'

        # Q
        Qi      =   inv(Znn * ZnnH)

        # 计算并写入结果
        preM[cbfs, neibfs] .=    view(ZnnH, cbfsInCnnei, :) * Qi

        # 更新进度条
        next!(pmeter)

    end # iCube
    # 恢复BLAS默认线程以防影响其他多线程函数
    BLAS.set_num_threads(nthds)
    # 保存预条件类型
    open(joinpath(SimulationParams.resultDir, "InputArgs.txt"), "a+")  do f
        write(f, "\npreT:\tSAI")
    end
    return SAIPrec{CT}(preM)
end


"""
计算稀疏近似逆 (Sparse Approximate Inverse (SAI)) 的函数
输入为近场阻抗矩阵CSC, 叶层信息 (也可以为非叶层, 但计算量更大) 
ZnearCSC::::SparseMatrixCSC{CT, Int}
cubes
该函数提供右预条件
"""
function sparseApproximateInversePr(ZnearCSC::ZnearT{CT}, cubes::AbstractVector) where { FT<:Real, CT<:Complex{FT}}
    # 将本函数内的BLAS设为单线程
    nthds = nthreads()
    BLAS.set_num_threads(1)
    # 首先预分配结果, 采用与 ZnearCSC 相同的稀疏模式
    preM    =   deepcopy(ZnearCSC)

    #TODO: 测试使用 CSR 格式矩阵做左预条件会不会快一点
    
    # 所有的盒子
    nCubes  =   length(cubes)

    # 进度条
    pmeter  =   Progress(nCubes, "Calculating SAI right preconditioner...")

    # 对所有盒子循环
    @threads for iCube in 1:nCubes

        # 本盒子与所有邻盒子id
        cube    =   cubes[iCube]
        ineiIDs =   cube.neighbors

        # 本盒子所有基函数的数量
        nNeibfs =   0
        for ineiID in ineiIDs
            nNeibfs    +=  length(cubes[ineiID].bfInterval)
        end
        # 本盒子所有基函数的 id
        neibfs  =   Vector{Int}(undef, nNeibfs)
        nNeibfs   =   0
        for ineiID in ineiIDs
            n   =   length(cubes[ineiID].bfInterval)
            neibfs[(nNeibfs + 1):(nNeibfs + n)]    .=  cubes[ineiID].bfInterval
            nNeibfs   +=  n
        end
        # 本盒子所有基函数的 id
        # cneibfs::Vector{Int}     =   vcat([collect(cubes[ineiID].bfInterval) for ineiID in ineiIDs]...)
        
        # 本盒子与所有邻盒子及其邻盒子的数量
        niNeiNeis=   0
        for ineiID in ineiIDs
            niNeiNeis    +=  length(cubes[ineiID].neighbors)
        end
        # 本盒子与所有邻盒子 id
        iNeisNeiIDs =   Vector{Int}(undef, niNeiNeis)
        niNeiNeis   =   0
        for ineiID in ineiIDs
            n   =   length(cubes[ineiID].neighbors)
            iNeisNeiIDs[(niNeiNeis + 1):(niNeiNeis + n)]    .=  cubes[ineiID].neighbors
            niNeiNeis   +=  n
        end
        unique!(sort!(iNeisNeiIDs))
        # iNeisNeiIDs::Vector{Int} =   unique(sort!(vcat([cubes[ineiID].neighbors for ineiID in ineiIDs]...)))

        # 本盒子所有基函数的数量
        nNeisNeibfs =   0
        for ineiID in iNeisNeiIDs
            nNeisNeibfs    +=  length(cubes[ineiID].bfInterval)
        end
        # 本盒子所有基函数的 id
        neisNeibfs  =   Vector{Int}(undef, nNeisNeibfs)
        nNeisNeibfs   =   0
        for ineiID in iNeisNeiIDs
            n   =   length(cubes[ineiID].bfInterval)
            neisNeibfs[(nNeisNeibfs + 1):(nNeisNeibfs + n)]    .=  cubes[ineiID].bfInterval
            nNeisNeibfs   +=  n
        end
        
        # 本盒子以及邻盒子的基函数 ids 数量
        # cnneibfs::Vector{Int}    =   vcat([collect(cubes[ineiID].bfInterval) for ineiID in iNeisNeiIDs]...)
        
        # 本盒子基函数
        cbfs        =   cube.bfInterval
        # 本盒子基函数起始点在 nneibfs 的位置
        cbfsInCnnei =   cbfs .+ (searchsortedfirst(neisNeibfs, cbfs.start) - cbfs.start)

        # 提取对应的阻抗矩阵
        Znn     =   zeros(CT, nNeisNeibfs, nNeibfs)
        Znn    .=   ZnearCSC[neisNeibfs, neibfs]
        ZnnH    =   adjoint(Znn)

        # bb      =   sparse(cbfsInCnnei, 1:length(cbfsInCnnei), one(CT), length(neisNeibfs), length(cbfsInCnei))

        # Q
        Qi      =   inv(ZnnH * Znn)

        # 计算并写入结果
        preM[neibfs, cbfs] .=    Qi * view(ZnnH, :, cbfsInCnnei)

        # 更新进度条
        next!(pmeter)

    end # iCube
    # 恢复BLAS默认线程以防影响其他多线程函数
    BLAS.set_num_threads(nthds)
    # 保存预条件类型
    open(MoM.SimulationParams.resultDir*"/InputArgs.txt", "a+")  do f
        write(f, "\npreT:\tSAI")
    end

    return SAIPrec{CT}(preM)
end

"""
计算稀疏近似逆 (Sparse Approximate Inverse (SAI)) 的函数
输入为近场阻抗矩阵CSC, 叶层信息 (也可以为非叶层, 但计算量更大) 
ZnearCSC::::SparseMatrixCSC{CT, Int}
level
该函数提供左预条件
"""
function sparseApproximateInversePl(ZnearCSC::ZnearT{CT}, level::LT) where { CT<:Complex, LT <: AbstractLevel}
    # 计算结果
    sparseApproximateInversePl(ZnearCSC, level.cubes)
end

"""
计算稀疏近似逆 (Sparse Approximate Inverse (SAI)) 的函数
输入为近场阻抗矩阵CSC, 叶层信息 (也可以为非叶层, 但计算量更大) 
ZnearCSC::::SparseMatrixCSC{CT, Int}
octree::Octree{FT}
该函数提供右预条件
"""
function sparseApproximateInversePr(ZnearCSC::ZnearT{CT}, level::LT) where { CT<:Complex, LT <: AbstractLevel}
    # 计算结果
    sparseApproximateInversePr(ZnearCSC, level.cubes)
end

"""
计算稀疏近似逆 (Sparse Approximate Inverse (SAI)) 的函数
输入为近场阻抗矩阵CSC, 叶层信息 (也可以为非叶层, 但计算量更大) 
ZnearCSC::::SparseMatrixCSC{CT, Int}
octree::Octree{FT}
该函数提供左预条件
"""
function sparseApproximateInversePl(ZnearCSC::ZnearT{CT}, octree::OctreeInfo{FT, LT}) where { FT<:Real, CT<:Complex{FT}, LT}
    # 层数
    nLevel::Int = max(keys(octree.levels)...)
    # 计算结果
    sparseApproximateInversePl(ZnearCSC, octree.levels[nLevel])
end

"""
计算稀疏近似逆 (Sparse Approximate Inverse (SAI)) 的函数
输入为近场阻抗矩阵CSC, 叶层信息 (也可以为非叶层, 但计算量更大) 
ZnearCSC::::SparseMatrixCSC{CT, Int}
octree::Octree{FT}
该函数提供右预条件
"""
function sparseApproximateInversePr(ZnearCSC::ZnearT{CT}, octree::OctreeInfo{FT, LT}) where { FT<:Real, CT<:Complex{FT}, LT}
    # 层数
    nLevel::Int = max(keys(octree.levels)...)
    # 计算结果
    sparseApproximateInversePr(ZnearCSC, octree.levels[nLevel])
end