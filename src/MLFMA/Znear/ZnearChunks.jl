abstract type ZNEARCHUNK{T} <: AbstractMatrix{T} end

"""
创建近场矩阵结构体，所包含的数据为所有盒子内的近场矩阵元，多线程版本
```
m::Int                  行数
n::Int                  列数
nChunks::Int            矩阵块儿数
chunks::Vector{ZnearChunksStruct{T}}    矩阵
lmul::Vector{T}         用于左乘其它矩阵、向量的临时数组，大小与列数相同
lmuld::Vector{T}        用于左乘其它矩阵、向量的临时分布式数组，大小与列数相同，默认不分配
rmul::Vector{T}         用于右乘其它矩阵、向量的临时数组，大小与行数相同
lmuld::Vector{T}        用于左乘其它矩阵、向量的临时分布式数组，大小与列数相同，默认不分配
```
"""
mutable struct ZnearChunksStruct{T<:Number} <: ZNEARCHUNK{T}
    m::Int
    n::Int
    nChunks::Int
    chunks::Vector{MatrixChunk{T}}
    lmul::Vector{T}
    lmuld::Vector{T}
    rmul::Vector{T}
    rmuld::Vector{T}
end

"""
    ZnearChunksStruct{T}(chunks; m, n) where {T<:Number}

ZnearChunksStruct 类的初始化函数。
"""
function ZnearChunksStruct{T}(chunks; m, n) where {T<:Number}

    nChunks =   length(chunks)
    lmul    =   zeros(T, n)
    lmuld   =   zeros(T, 0)
    rmul    =   zeros(T, m)
    rmuld   =   zeros(T, 0)

    ZnearChunksStruct{T}(m, n, nChunks, chunks, lmul, lmuld, rmul, rmuld)
end


import Base:*, size, show, display, eltype, length, getindex, setindex!

size(Z::T) where{T<:ZNEARCHUNK}  = (Z.m, Z.n)
size(Z::T, d::Integer) where{T<:ZNEARCHUNK}  = size(Z)[d]
# show(Z::T) where{T<:ZNEARCHUNK}  = show(Z.mat)
show(io::IO, unused::MIME{Symbol("text/plain")}, Z::T) where{T<:ZNEARCHUNK} = begin
    println("A $(Z.m) × $(Z.n) matrix made of chunks.")
end

display(Z::T) where{T<:ZNEARCHUNK} = begin
    println("A $(Z.m) × $(Z.n) matrix made of chunks.")
end

eltype(Z::T) where{T<:ZNEARCHUNK} = eltype(Z.lmul)
length(Z::T) where{T<:ZNEARCHUNK} = Z.m*Z.n

"""
重载 getindex 函数
"""
function getindex(Z::T, i1::Int, i2::Int) where{T<:ZNEARCHUNK}
    
    chunki  =   searchsortedfirst(Z.chunks, i1; lt = (ck, i1) -> first(ck.rowIndices) <= i1)

    Zchunk  =   Z.chunks[chunki]
    mptr =   searchsorted(Zchunk.rowIndices, i1)
    nptr =   searchsorted(Zchunk.colIndices, i2)
    (isempty(mptr) || isempty(nptr)) && return zero(eltype(Zchunk.mat))
    return Zchunk.mat[mptr[1], nptr[1]]

end

"""
    setindex!(Z::T, x, i1::Int, i2::Int) where{T<:ZNEARCHUNK}

重载 setindex! 函数
"""
function setindex!(Z::T, x, i1::Int, i2::Int) where{T<:ZNEARCHUNK}

    chunki  =   searchsortedfirst(Z.chunks, i1; lt = (ck, i1) -> first(ck.rowIndices) <= i1)

    Zchunk  =   Z.chunks[chunki]
    mptr =   searchsorted(Zchunk.rowIndices, i1)
    nptr =   searchsorted(Zchunk.colIndices, i2)
    (isempty(mptr) || isempty(nptr)) && throw("目标索引不在本区间内")
    Zchunk.mat[mptr[1], nptr[1]] = x
    return nothing
end


"""
    initialZchunksMulV!(Z::T) where{T<:ZnearChunksStruct}

初始化 阻抗矩阵 左乘 向量 乘积的 分布式数组。
"""
function initialZchunksMulV!(Z::T) where{T<:ZnearChunksStruct}

    Z.rmuld = Z.rmul

    nothing
end

"""
    initialVMulZchunks!(Z::T) where{T<:ZnearChunksStruct}

初始化 阻抗矩阵 右乘 向量 乘积的 分布式数组
"""
function initialVMulZchunks!(Z::T) where{T<:ZnearChunksStruct}

    Z.lmuld = Z.lmul
    nothing
end


"""
    Base.:*(Z::T, x::AbstractVector) where{T<:ZnearChunksStruct}

实现左乘其它向量
"""
function Base.:*(Z::T, x::AbstractVector) where{T<:ZnearChunksStruct}
    # initialZchunksMulV!(Z)
    nthds   =   nthreads()
    BLAS.set_num_threads(1)
    @inbounds @threads  for ii in 1:Z.nChunks
        copyto!(view(Z.rmul, Z.chunks[ii].rowIndices), Z.chunks[ii] * x)
    end
    BLAS.set_num_threads(nthds)

    return copy(Z.rmul)
end

function LinearAlgebra.mul!(y, Z::T, x::AbstractVector) where{T<:ZnearChunksStruct}
    # initialZchunksMulV!(Z)
    nthds   =   nthreads()
    BLAS.set_num_threads(1)
    @inbounds @threads  for ii in 1:Z.nChunks
        mul!(view(Z.rmul, Z.chunks[ii].rowIndices), Z.chunks[ii], x)
    end
    BLAS.set_num_threads(nthds)

    copyto!(y, Z.rmul)

    return y
end


"""
    get_chunks_minmax_col(matchunks)


"""
function get_chunks_minmax_col(matchunks)::UnitRange{Int64}

    mincol = minimum(matchunk -> first(matchunk.colIndices), matchunks)
    maxcol = maximum(matchunk -> last(matchunk.colIndices), matchunks)

    return mincol:maxcol
    
end

"""
    Base.:*(Z::ZNEARCHUNK{T}, mat::AbstractMatrix) where{T<:Number}

实现左乘其它矩阵
"""
function Base.:*(Z::ZNEARCHUNK{T}, mat::AbstractMatrix) where{T<:Number}
    initialZchunksMulV!(Z)
    nre =   size(mat, 2)
    re  =   zeros(T, Z.m, size(mat, 2))
    @inbounds for i in 1:nre
        re[:, i] .= Z * view(mat, :, i)
    end

    return re
end

