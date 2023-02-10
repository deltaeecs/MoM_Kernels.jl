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

    # 所有的盒子
    nCubes  =   length(cubes)

    # 进度条
    pmeter  =   Progress(nCubes, "Pₗ")

    # Znn ZnnH 按线程预分配内存
    Znnts       =   [zeros(CT, 1) for _ in 1:nthds]
    ZnnHZnnts   =   [zeros(CT, 1) for _ in 1:nthds]

    # 对所有盒子循环
    @threads for iCube in 1:nCubes

        # 本盒子与所有邻盒子id
        cube    =   cubes[iCube]
        ineiIDs =   cube.neighbors

        # 本盒子与所有邻盒子基函数的 id
        neibfs  = reduce(vcat, map(i -> cubes[i].bfInterval, ineiIDs))
        # 本盒子所有邻盒子基函数的数量
        nNeibfs = length(neibfs)

        # 本盒子与所有邻及其邻盒子 id
        iNeisNeiIDs  = reduce(vcat, map(i -> cubes[i].neighbors, ineiIDs))
        unique!(sort!(iNeisNeiIDs))

        # 本盒子与所有邻盒子及其邻盒子基函数的 id
        neisNeibfs  = reduce(vcat, map(i -> cubes[i].bfInterval, iNeisNeiIDs))
        # 本盒子与所有邻盒子及其邻盒子基函数的数量
        nNeisNeibfs = length(neisNeibfs)

        # 提取对应的阻抗矩阵
        # Znn     =   zeros(CT, nNeibfs, nNeisNeibfs)
        # Znn 保存在预分配内存里
        Znnt    =   Znnts[threadid()]
        length(Znnt) < nNeibfs*nNeisNeibfs && resize!(Znnt, nNeibfs*nNeisNeibfs)
        Znn     =   reshape(view(Znnt, 1:nNeibfs*nNeisNeibfs), nNeibfs, nNeisNeibfs)

        for j in 1:length(neisNeibfs), i in 1:length(neibfs)
            Znn[i, j]  =   ZnearCSC[neibfs[i], neisNeibfs[j]]
        end
        ZnnH    =   Znn'

        # Znn*ZnnH 也保存在预分配内存里
        ZnnHZnnt = ZnnHZnnts[threadid()]
        length(ZnnHZnnt) < nNeibfs*nNeibfs && resize!(ZnnHZnnt, nNeibfs*nNeibfs)
        ZnnHZnn     =   reshape(view(ZnnHZnnt, 1:nNeibfs*nNeibfs), nNeibfs, nNeibfs)

        # Qi      =   inv(Znn * ZnnH)
        Qi      =   inv(mul!(ZnnHZnn, Znn, ZnnH))

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

    # 所有的盒子
    nCubes  =   length(cubes)

    # 进度条
    pmeter  =   Progress(nCubes, "Calculating SAI right preconditioner...")

    # Znn ZnnH 按线程预分配内存
    Znnts       =   [zeros(CT, 1) for _ in 1:nthds]
    ZnnZnnHts   =   [zeros(CT, 1) for _ in 1:nthds]

    # 对所有盒子循环
    @threads for iCube in 1:nCubes

        # 本盒子与所有邻盒子id
        cube    =   cubes[iCube]
        ineiIDs =   cube.neighbors

        # 本盒子与所有邻盒子基函数的 id
        neibfs  = reduce(vcat, map(i -> cubes[i].bfInterval, ineiIDs))
        # 本盒子所有邻盒子基函数的数量
        nNeibfs = length(neibfs)

        # 本盒子与所有邻及其邻盒子 id
        iNeisNeiIDs  = reduce(vcat, map(i -> cubes[i].neighbors, ineiIDs))
        unique!(sort!(iNeisNeiIDs))

        # 本盒子与所有邻盒子及其邻盒子基函数的 id
        neisNeibfs  = reduce(vcat, map(i -> cubes[i].bfInterval, iNeisNeiIDs))
        # 本盒子与所有邻盒子及其邻盒子基函数的数量
        nNeisNeibfs = length(neisNeibfs)

        # 提取对应的阻抗矩阵
        # Znn 保存在预分配内存里
        Znnt    =   Znnts[threadid()]
        length(Znnt) < nNeibfs*nNeisNeibfs && resize!(Znnt, nNeibfs*nNeisNeibfs)
        Znn     =   reshape(view(Znnt, 1:nNeibfs*nNeisNeibfs), nNeisNeibfs, nNeibfs)
        # Znn    .=   ZnearCSC[neisNeibfs, neibfs]
        for j in 1:length(neibfs), i in 1:length(neisNeibfs)
            Znn[i, j]  =   ZnearCSC[neisNeibfs[i], neibfs[j]]
        end
        ZnnH    =   adjoint(Znn)

        # Znn*ZnnH 也保存在预分配内存里
        ZnnZnnHt = ZnnZnnHts[threadid()]
        length(ZnnZnnHt) < nNeibfs*nNeibfs && resize!(ZnnZnnHt, nNeibfs*nNeibfs)
        ZnnZnnH     =   reshape(view(ZnnZnnHt, 1:nNeibfs*nNeibfs), nNeibfs, nNeibfs)

        # Qi      =   inv(ZnnH * Znn)
        Qi      =   inv(mul!(ZnnZnnH, ZnnH, Znn))

        # 计算并写入结果
        preM[neibfs, cbfs] .=    Qi * view(ZnnH, :, cbfsInCnnei)

        # 更新进度条
        next!(pmeter)

    end # iCube
    # 恢复BLAS默认线程以防影响其他多线程函数
    BLAS.set_num_threads(nthds)
    # 保存预条件类型
    open(SimulationParams.resultDir*"/InputArgs.txt", "a+")  do f
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