# 本文件提供计算稀疏近似逆预条件的函数
"""
    SAIChunkPrec{T} <: AbstractMatrix{T}

分块系数近似逆的结构体。
"""
struct SAIChunkPrec{T} <: AbstractMatrix{T}
    mat::ZNEARCHUNK{T}
end

Base.show(io::IO, unused::MIME{Symbol("text/plain")}, M::SAIChunkPrec{T}) where {T} = show(io, unused, M.mat)
Base.display(M::SAIChunkPrec{T}) where {T} = display(M.mat)
Base.size(M::SAIChunkPrec{T}) where {T} = size(M.mat)

# nnz(M::SAIChunkPrec{T}) where {T} = nnz(M.mat)

@doc """
    ldiv!(M::SAIChunkPrec{T}, x::AbstractVector) where {T}
    ldiv!(y::AbstractVector, M::SAIChunkPrec{T}, x::AbstractVector) where {T}

实现 `x .= M * x` 或 `y .= M * x`。
"""
function LinearAlgebra.ldiv!(M::SAIChunkPrec{T}, x::AbstractVector) where {T}
    mul!(x, M.mat, x)
end
function LinearAlgebra.ldiv!(y::AbstractVector, M::SAIChunkPrec{T}, x::AbstractVector) where {T}
    mul!(y, M.mat, x)
end

Base.eltype(::SAIChunkPrec{T}) where {T} = T
Base.:\(M::SAIChunkPrec{T}, x::AbstractVector) where {T}= ldiv!(copy(x), M, x)

"""
    sparseApproximateInversePl(ZnearChunks::ZnearChunksStruct{CT}, level; nbf = 0) where {FT<:Real, CT<:Complex{FT}}

根据块状近场阻抗矩阵 `ZnearChunks` 和计算阻抗矩阵层的盒子信息 `cubes` 计算左稀疏近似逆 (Sparse Approximate Inverse (SAI)) 。
"""
function sparseApproximateInversePl(ZnearChunks::ZnearChunksStruct{CT}, level; nbf = 0) where {FT<:Real, CT<:Complex{FT}}

    # 将本函数内的BLAS设为单线程
    nthds = nthreads()
    BLAS.set_num_threads(1)
    # 首先预分配结果, 采用与 Znear 相同的稀疏模式
    preM    =   deepcopy(ZnearChunks)
    
    cubes   =   level.cubes
    # 所有的盒子
    nCubes  =   length(cubes)

    # 进度条
    pmeter  =   Progress(nCubes; desc = "Pₗ (T) ...")
    
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

        Znn     =   zeros(CT, nNeibfs, nNeisNeibfs)
        # 提取对应的阻抗矩阵
        for cubeNeiNei in ineiIDs
            # cubeNeiNei = 3
            ZChunkNeiNei =  ZnearChunks.chunks[cubeNeiNei]
            rowIndices   =  ZChunkNeiNei.rowIndices
            colIndices   =  ZChunkNeiNei.colIndices
            for (i, m) in enumerate(rowIndices)
                miZnn   =   searchsortedfirst(neibfs, m)
                for (j, n) in enumerate(colIndices)
                    njZnn   =   searchsortedfirst(neisNeibfs, n)
                    Znn[miZnn, njZnn] = ZChunkNeiNei.mat[i, j]
                end
            end
        end

        ZnnH    =   Znn'
        
        # cbfsInCnnei =   cbfs .+ (searchsortedfirst(neibfs, cbfs.start) - cbfs.start)

        # Znn     =   zeros(CT, nNeibfs, nNeibfs)
        # for j in 1:length(neibfs), i in 1:length(neibfs)
        #     Znn[i, j]  =   ZnearChunks[neibfs[i], neibfs[j]]
        # end
        # ZnnH    =   Znn'

        # Q
        Qi      =   inv(Znn * ZnnH)

        # 计算并写入结果
        preM.chunks[iCube].mat .=    view(ZnnH, cbfsInCnnei, :) * Qi

        # 更新进度条
        next!(pmeter)

    end # iCube
    # 恢复BLAS默认线程以防影响其他多线程函数
    BLAS.set_num_threads(nthds)
    # 保存预条件类型
    open(joinpath(SimulationParams.resultDir, "InputArgs.txt"), "a+")  do f
        @printf f "%20s\n" "预条件"
        @printf f "%-20s %13s\n" "类型" "SAI" 
    end
    return SAIChunkPrec{CT}(preM)
end