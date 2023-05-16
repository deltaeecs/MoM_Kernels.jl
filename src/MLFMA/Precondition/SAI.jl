# 本文件提供计算稀疏近似逆预条件的函数

"""
    SAIPrec{T} <: AbstractMatrix{T} 

封装系数近似逆矩阵的类型。
"""
struct SAIPrec{T} <: AbstractMatrix{T}
    mat::AbstractMatrix{T}
end

Base.show(io::IO, unused::MIME{Symbol("text/plain")}, M::SAIPrec{T}) where {T} = show(io, unused, M.mat)
Base.display(M::SAIPrec{T}) where {T} = display(M.mat)
Base.size(M::SAIPrec{T}) where {T} = size(M.mat)
SparseArrays.nnz(M::SAIPrec{T}) where {T} = nnz(M.mat)


@doc """
    ldiv!(M::SAIPrec{T}, x::AbstractVector) where {T}
    ldiv!(y::AbstractVector, M::SAIPrec{T}, x::AbstractVector) where {T}

实现 `x .= M * x` 或 `y .= M * x`。
"""
function LinearAlgebra.ldiv!(M::SAIPrec{T}, x::AbstractVector) where {T}
    x .= M.mat * x
end
function LinearAlgebra.ldiv!(y::AbstractVector, M::SAIPrec{T}, x::AbstractVector) where {T}
    mul!(y, M.mat, x)
end

Base.eltype(::SAIPrec{T}) where {T} = T
Base.:\(M::SAIPrec{T}, x::AbstractVector) where {T}= ldiv!(copy(x), M, x)

"""
    sparseApproximateInversePl(Znear::ZnearT{CT}, cubes::AbstractVector) where {FT<:Real, CT<:Complex{FT}}

根据近场阻抗矩阵 `Znear` 和计算阻抗矩阵层的盒子信息 `cubes` 计算左稀疏近似逆 (Sparse Approximate Inverse (SAI)) 。
"""
function sparseApproximateInversePl(Znear::ZnearT{CT}, cubes::AbstractVector) where {FT<:Real, CT<:Complex{FT}}
    @clock "计算SAI" begin
        # 将本函数内的BLAS设为单线程
        nthds = nthreads()
        BLAS.set_num_threads(1)
        # 首先预分配结果, 采用与 Znear 相同的稀疏模式
        preM    =   deepcopy(Znear)

        # 所有的盒子
        nCubes  =   length(cubes)

        # 进度条
        pmeter  =   Progress(nCubes, "Pₗ")

        # Znn ZnnH 按线程预分配内存
        Znnts       =   [zeros(CT, 1) for _ in 1:nthds]
        ZnnHZnnts   =   [zeros(CT, 1) for _ in 1:nthds]
        PHts        =   [zeros(CT, 1) for _ in 1:nthds]

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

            # 本盒子基函数
            cbfs        =   cube.bfInterval
            # 本盒子基函数起始点在 nneibfs 的位置
            cbfsInCnnei =   cbfs .+ (searchsortedfirst(neisNeibfs, cbfs.start) - cbfs.start)
            
            # 提取对应的阻抗矩阵
            # Znn     =   zeros(CT, nNeibfs, nNeisNeibfs)
            # Znn 保存在预分配内存里
            Znnt    =   Znnts[threadid()]
            length(Znnt) < nNeibfs*nNeisNeibfs && resize!(Znnt, nNeibfs*nNeisNeibfs)
            fill!(Znnt, 0)
            Znn     =   reshape(view(Znnt, 1:nNeibfs*nNeisNeibfs), nNeibfs, nNeisNeibfs)

            for j in 1:length(neisNeibfs), i in 1:length(neibfs)
                Znn[i, j]  =   Znear[neibfs[i], neisNeibfs[j]]
            end
            ZnnH    =   Znn'

            # Znn*ZnnH 也保存在预分配内存里
            ZnnHZnnt = ZnnHZnnts[threadid()]
            length(ZnnHZnnt) < nNeibfs*nNeibfs && resize!(ZnnHZnnt, nNeibfs*nNeibfs)
            fill!(ZnnHZnnt, 0)
            ZnnHZnn     =   reshape(view(ZnnHZnnt, 1:nNeibfs*nNeibfs), nNeibfs, nNeibfs)

            PHt = PHts[threadid()]
            length(PHt) < length(cbfs)*nNeibfs && resize!(PHt, length(cbfs)*nNeibfs)
            fill!(PHt, 0)
            PH  =   reshape(view(PHt, 1:length(cbfs)*nNeibfs), nNeibfs, length(cbfs))

            # Qi      =   inv(Znn * ZnnH)
            # 先计算出 ZnnHZnn， 然后对其进行 QR分解
            mul!(ZnnHZnn, Znn, ZnnH)
            ZLU      =   lu!(ZnnHZnn)

            # 将结果先写入 PHt
            ldiv!(PH, ZLU, view(Znn, :, cbfsInCnnei))

            # 写入结果
            preM[cbfs, neibfs] .=   PH'

            # 更新进度条
            next!(pmeter)

        end # iCube
        # 恢复BLAS默认线程以防影响其他多线程函数
        BLAS.set_num_threads(nthds)
    end
    # 保存预条件类型
    open(joinpath(SimulationParams.resultDir, "InputArgs.txt"), "a+")  do f
        @printf f "%20s\n" "预条件"
        @printf f "%-20s %13s\n" "类型" "SAI" 
    end
    return SAIPrec{CT}(preM)
end


"""
    sparseApproximateInversePr(Znear::ZnearT{CT}, cubes::AbstractVector) where {FT<:Real, CT<:Complex{FT}}

根据近场阻抗矩阵 `Znear` 和计算阻抗矩阵层的盒子信息 `cubes` 计算右稀疏近似逆 (Sparse Approximate Inverse (SAI)) 。
"""
function sparseApproximateInversePr(Znear::ZnearT{CT}, cubes::AbstractVector) where { FT<:Real, CT<:Complex{FT}}
    @clock "计算SAI" begin
        # 将本函数内的BLAS设为单线程
        nthds = nthreads()
        BLAS.set_num_threads(1)
        # 首先预分配结果, 采用与 Znear 相同的稀疏模式
        preM    =   deepcopy(Znear)

        # 所有的盒子
        nCubes  =   length(cubes)

        # 进度条
        pmeter  =   Progress(nCubes, "Calculating SAI right preconditioner...")

        # Znn ZnnH 按线程预分配内存
        Znnts       =   [zeros(CT, 1) for _ in 1:nthds]
        ZnnZnnHts   =   [zeros(CT, 1) for _ in 1:nthds]
        Pts         =   [zeros(CT, 1) for _ in 1:nthds]

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

            # 本盒子基函数
            cbfs        =   cube.bfInterval
            # 本盒子基函数起始点在 nneibfs 的位置
            cbfsInCnnei =   cbfs .+ (searchsortedfirst(neisNeibfs, cbfs.start) - cbfs.start)
            
            # 提取对应的阻抗矩阵
            # Znn 保存在预分配内存里
            Znnt    =   Znnts[threadid()]
            length(Znnt) < nNeibfs*nNeisNeibfs && resize!(Znnt, nNeibfs*nNeisNeibfs)
            fill!(Znnt, 0)
            Znn     =   reshape(view(Znnt, 1:nNeibfs*nNeisNeibfs), nNeisNeibfs, nNeibfs)
            # Znn    .=   Znear[neisNeibfs, neibfs]
            for j in 1:length(neibfs), i in 1:length(neisNeibfs)
                Znn[i, j]  =   Znear[neisNeibfs[i], neibfs[j]]
            end
            ZnnH    =   adjoint(Znn)

            Pt = Pts[threadid()]
            length(Pt) < length(cbfs)*nNeibfs && resize!(Pt, length(cbfs)*nNeibfs)
            fill!(Pt, 0)
            P  =   reshape(view(Pt, 1:length(cbfs)*nNeibfs), nNeibfs, length(cbfs))

            # Znn*ZnnH 也保存在预分配内存里
            ZnnZnnHt = ZnnZnnHts[threadid()]
            length(ZnnZnnHt) < nNeibfs*nNeibfs && resize!(ZnnZnnHt, nNeibfs*nNeibfs)
            fill!(ZnnZnnHt, 0)
            ZnnZnnH     =   reshape(view(ZnnZnnHt, 1:nNeibfs*nNeibfs), nNeibfs, nNeibfs)

            # 先计算出 ZnnHZnn， 然后对其进行 QR分解
            mul!(ZnnZnnH, ZnnH, Znn)
            ZLU      =   lu!(ZnnZnnH)

            # 将结果先写入 P
            ldiv!(P, ZLU, view(ZnnH, :, cbfsInCnnei))

            # 写入结果
            preM[neibfs, cbfs] .=   P

            # 更新进度条
            next!(pmeter)

        end # iCube
        # 恢复BLAS默认线程以防影响其他多线程函数
        BLAS.set_num_threads(nthds)
    end
    # 保存预条件类型
    open(joinpath(SimulationParams.resultDir, "InputArgs.txt"), "a+")  do f
        @printf f "%20s\n" "预条件"
        @printf f "%-20s %13s\n" "类型" "SAI" 
    end

    return SAIPrec{CT}(preM)
end

"""
    sparseApproximateInversePl(Znear::ZnearT{CT}, level::LT) where { CT<:Complex, LT <: AbstractLevel}

根据近场阻抗矩阵 `Znear` 和计算阻抗矩阵层信息 `level` 计算左稀疏近似逆 (Sparse Approximate Inverse (SAI)) 。
"""
function sparseApproximateInversePl(Znear::ZnearT{CT}, level::LT) where { CT<:Complex, LT <: AbstractLevel}
    # 计算结果
    sparseApproximateInversePl(Znear, level.cubes)
end

"""
    sparseApproximateInversePr(Znear::ZnearT{CT}, level::LT) where { CT<:Complex, LT <: AbstractLevel}

根据近场阻抗矩阵 `Znear` 和计算阻抗矩阵层信息 `level` 计算右稀疏近似逆 (Sparse Approximate Inverse (SAI)) 。
"""
function sparseApproximateInversePr(Znear::ZnearT{CT}, level::LT) where { CT<:Complex, LT <: AbstractLevel}
    # 计算结果
    sparseApproximateInversePr(Znear, level.cubes)
end

"""
    sparseApproximateInversePl(Znear::ZnearT{CT}, octree::OctreeInfo{FT, LT}) where { FT<:Real, CT<:Complex{FT}, LT}
    
根据近场阻抗矩阵 `Znear` 和八叉树 `octree` 叶层计算左稀疏近似逆 (Sparse Approximate Inverse (SAI)) 。
"""
function sparseApproximateInversePl(Znear::ZnearT{CT}, octree::OctreeInfo{FT, LT}) where { FT<:Real, CT<:Complex{FT}, LT}
    # 层数
    nLevel::Int = max(keys(octree.levels)...)
    # 计算结果
    sparseApproximateInversePl(Znear, octree.levels[nLevel])
end

"""
    sparseApproximateInversePr(Znear::ZnearT{CT}, octree::OctreeInfo{FT, LT}) where { FT<:Real, CT<:Complex{FT}, LT}

根据近场阻抗矩阵 `Znear` 和八叉树 `octree` 叶层计算右稀疏近似逆 (Sparse Approximate Inverse (SAI)) 。
"""
function sparseApproximateInversePr(Znear::ZnearT{CT}, octree::OctreeInfo{FT, LT}) where { FT<:Real, CT<:Complex{FT}, LT}
    # 层数
    nLevel::Int = max(keys(octree.levels)...)
    # 计算结果
    sparseApproximateInversePr(Znear, octree.levels[nLevel])
end