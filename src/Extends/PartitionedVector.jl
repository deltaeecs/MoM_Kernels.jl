
"""
    PartitionedVector{T} <: AbstractVector{T}

用于保存向量块的类，同时在块内保存一些其他块的数据。

```
size::Int                           原始 Vector 的大小
data::OffsetVector{T, Vector{T}}    本地保存的数据
indices::UnitRange{Int}             本地保存数据的索引区间
ghostdata::SparseVector{T, Int}     用到的其它数据
ghostindices::Vector{Int}           用到的其它数据的索引区间
```
"""
struct PartitionedVector{T} <: AbstractVector{T}
    size::Int
    data::OffsetVector{T, Vector{T}}
    indices::UnitRange{Int}
    ghostdata::SparseVector{T, Int}
    ghostindices::Vector{Int}
end

Base.size(A::T) where{T<:PartitionedVector}  = (A.size, )
Base.eltype(::PartitionedVector{T}) where{T} = T
Base.length(A::T) where{T<:PartitionedVector} = A.size

"""
    Base.getindex(A::PartitionedVector, i::I) where {I<:Integer}

重载 `PartitionedVector` 类型的 `getindex` 函数。
"""
function Base.getindex(A::PartitionedVector, i::I) where {I<:Integer}
    if i in A.indices
        return getindex(A.data, i)
    elseif !isempty(searchsorted(A.ghostindices, i))
        return getindex(A.ghostdata, i)
    else
        @warn "$i is out of indices of Array"
		return nothing
    end
end