"""
创建近场矩阵块结构体，所包含的数据为某一盒子内的近场矩阵元
```
m::Int，行数
n::Int，列数
mat::Matrix{T}，矩阵
rowIndices::AbstractVector{Int}，行索引
colIndices::AbstractVector{Int}，列索引
lmul::AbstractVector{T}，用于左乘其它矩阵、向量的临时数组，大小与列数相同
rmul::AbstractVector{T}，用于右乘其它矩阵、向量的临时数组，大小与行数相同
```
"""
struct MatrixChunk{T<:Number} <:AbstractMatrix{T}
    m::Int
    n::Int
    mat::Matrix{T}
    rowIndices::UnitRange{Int}
    colIndices::Vector{Int}
    lmul::Vector{T}
    rmul::Vector{T}
end


function MatrixChunk{T}(rowIndices, colIndices) where {T<:Number}
    m   =   length(rowIndices)
    n   =   length(colIndices)
    mat =   zeros(T, m, n)
    lmul=   zeros(T, n)
    rmul=   zeros(T, m)

    MatrixChunk{T}(m, n, mat, rowIndices, colIndices, lmul, rmul)
end

import Base:*, size, show, display, eltype, length, getindex, setindex!

size(Z::T) where{T<:MatrixChunk}  = size(Z.mat)
size(Z::T, d::Integer) where{T<:MatrixChunk}  = size(Z.mat, d)
# show(Z::T) where{T<:MatrixChunk}  = show(Z.mat)
show(io::IO, unused::MIME{Symbol("text/plain")}, Z::T) where{T<:MatrixChunk} = begin
    println("A $(Z.m) × $(Z.n) matrix chunk.")
end
display(Z::T) where{T<:MatrixChunk} = begin
    println("A $(Z.m) × $(Z.n) matrix chunk.")
end
eltype(Z::T) where{T<:MatrixChunk} = eltype(Z.mat)
length(Z::T) where{T<:MatrixChunk} = Z.m*Z.n

"""
重载 getindex 函数
"""
function getindex(Z::T, i1::Int, i2::Int) where{T<:MatrixChunk}
    
    mptr =   searchsorted(Z.rowIndices, i1)
    nptr =   searchsorted(Z.colIndices, i2)
    (isempty(mptr) || isempty(nptr)) && return zero(eltype(Z.mat))
    return Z.mat[mptr[1], nptr[1]]

end

"""
重载 setindex! 函数
"""
function setindex!(Z::T, x, i1::Int, i2::Int) where{T<:MatrixChunk}

    mptr =   searchsorted(Z.rowIndices, i1)
    nptr =   searchsorted(Z.colIndices, i2)

    (isempty(mptr) || isempty(nptr)) && throw("目标索引不在本区间内")
    Z.mat[mptr[1], nptr[1]] = x
    return nothing
end


"""
实现左乘其它向量，默认矩阵块较小，不在本阶段并行
"""
function Base.:*(Z::T, x::AbstractVector) where{T<:MatrixChunk}
    
    @inbounds for i in 1:Z.n
        Z.lmul[i] = x[Z.colIndices[i]]
    end

    mul!(Z.rmul, Z.mat, Z.lmul)

    return Z.rmul
end

function LinearAlgebra.mul!(y, Z::T, x::AbstractVector) where{T<:MatrixChunk}
    
    @inbounds for i in 1:Z.n
        Z.lmul[i] = x[Z.colIndices[i]]
    end

    mul!(y, Z.mat, Z.lmul)

    return y
end

"""
实现左乘其它矩阵，默认矩阵块较小，不在本阶段并行
"""
function Base.:*(Z::T, mat::AbstractMatrix) where{T<:MatrixChunk}

    nre =   size(mat, 2)
    re  =   zeros(T, Z.m, size(mat, 2))
    @inbounds for i in 1:nre
        re[:, i] .= Z * view(mat, :, i)
    end

    return re
end

"""
实现右乘其它向量，默认矩阵块较小，不在本阶段并行
"""
function Base.:*(x::AbstractVector, Z::T) where{T<:MatrixChunk}
    
    @inbounds for i in 1:Z.m
        Z.rmul[i] = x[Z.rowIndices[i]]
    end

    copyto!(Z.lmul, Z.rmul * Z.mat)

    return Z.lmul
end

"""
实现右乘其它矩阵，默认矩阵块较小，不在本阶段并行
"""
function Base.:*(mat::AbstractMatrix, Z::T) where{T<:MatrixChunk}

    mre =   size(mat, 1)
    re  =   zeros(T, mre, Z.n)
    @inbounds for i in 1:mre
        re[i, :] .= view(mat, i, :) * Z
    end

    return re
end

"""
存储分布式矩阵向量乘积结果的数组
"""
struct mmulvStruct{T<:Number} <: AbstractVector{T}
    indices::AbstractVector{Int}
    vals::Vector{T}
end

eltype(v::mmulvStruct) = eltype(v.vals)
size(v::mmulvStruct) = size(v.vals)
length(v::mmulvStruct) = length(v.vals)

"""
重载 getindex 函数
"""
function getindex(v::T, i1::Int) where{T<:mmulvStruct}
    
    mptr =   searchsorted(v.indices, i1)
    isempty(mptr) && return zero(eltype(v))
    return v.vals[mptr[1]]

end

"""
重载 setindex! 函数
"""
function setindex!(v::T, x, i1::Int) where{T<:mmulvStruct}

    mptr =   searchsorted(v.indices, i1)

    isempty(mptr) && throw("目标索引不在本区间内")
    v.vals[mptr[1]] = x
    return nothing
end

"""
根据八叉树盒子信息初始化 cube 对应的近场矩阵元块儿
"""
function initialmmulvStruct(ZnearChunk)
    mmulvStruct{eltype(ZnearChunk)}(ZnearChunk.rowIndices, ZnearChunk.rmul)
end

"""
根据八叉树盒子信息初始化 cube 对应的近场矩阵元块儿
"""
function initialvmulmStruct(ZnearChunk)
    mmulvStruct{eltype(ZnearChunk)}(ZnearChunk.rowIndices, ZnearChunk.lmul)
end