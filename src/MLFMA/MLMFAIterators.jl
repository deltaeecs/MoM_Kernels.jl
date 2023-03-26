## 本文件的函数用于构建迭代算子
include("AggregateOnLevel.jl")

"""
保存 MLFMA 相关信息的结构体
"""
mutable struct MLMFAIterator{T<:Number, VT<:AbstractVector}
    octree          ::OctreeInfo
    Znear           ::AbstractMatrix{T}
    vsCellsInfo     ::CellT where {CellT<:AbstractVector}
    bfsInfo         ::BFT where {BFT<:AbstractVector}
    aggSBF          ::AbstractArray{T}
    disaggSBF       ::AbstractArray{T}
    ZI              ::VT
end

# 各种可作为线性算子的集合
LinearMapType{T}   =   Union{AbstractMatrix{T}, MLMFAIterator{T, VT}, LinearMap{T}} where {T<:Number, VT}

Base.eltype(::MLMFAIterator{T, VT}) where {T, VT} = T
Base.size(opt::MLMFAIterator) = size(opt.Znear)
Base.size(opt::MLMFAIterator, ind...) = size(opt.Znear, ind...)

function Base.getproperty(obj::MLMFAIterator, sym::Symbol)
    if sym === :leafLevel
        return obj.octree.levels[obj.nLevels]
    elseif sym === :nLevels
        return obj.octree.nLevels
    elseif sym === :levels
        return obj.octree.levels
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

# 多线程迭代函数
include("IterateOnOctree.jl")


"""
实现矩阵向量乘积，并封装为线性算子
"""
function MLMFAIterator(ZnearCSC, octree::OctreeInfo{FT, LT}, 
    vsCellsInfo::Vector, bfsInfo::Vector) where {FT<:Real, LT<:LevelInfo}
    
    # 叶层ID
    nLevels     =   octree.nLevels
    # 层信息
    levels      =   octree.levels
    # 叶层
    leafLevel   =   octree.levels[nLevels]
    
    # 基函数数量
    nbf         =   getNUnknown(bfsInfo)
    
    # 预先计算叶层聚合项，内存占用为该层采样点数 nPoles × Nbf， 因此仅在内存充足时使用
    @clock "叶层聚合项计算" begin
        aggSBF, disaggSBF   =   getAggSBFOnLevel(leafLevel, vsCellsInfo, bfsInfo)
    end

    # 给各层的聚合项、解聚项预分配内存
    memoryAllocationOnLevels!(nLevels, levels)

    # 给矩阵向量乘积预分配内存
    ZI  =   zeros(Complex{FT}, nbf)

    Zopt =  MLMFAIterator{Complex{FT}, Vector}(octree, ZnearCSC, vsCellsInfo, bfsInfo, aggSBF, disaggSBF, ZI)

    record_memorys(Zopt)

    return Zopt

end

"""
    LinearAlgebra.mul!(y, Zopt::MLMFAIterator, x)

    重载以实现矩阵向量乘积计算
TBW
"""
function LinearAlgebra.mul!(y, Zopt::MLMFAIterator, x)
    # near
    fill!(y, 0)
    mul!(y, Zopt.Znear, x)
    # far
    calZfarI!(Zopt, x)
    # cumsum results
    axpy!(1, Zopt.ZI, y)

    return y

end

"""
    LinearAlgebra.mul!(y, Zopt::MLMFAIterator, x)

    重载以实现矩阵向量乘积计算
TBW
"""
function LinearAlgebra.mul!(y::AbstractVector, Zopt::MLMFAIterator, x::AbstractVector, α::Number, β::Number)
    # near
    mul!(y, Zopt.Znear, x, α, β)
    # far
    calZfarI!(Zopt, x)
    # cumsum results
    axpy!(α, Zopt.ZI, y)

    return y

end

Base.:*(opt::MLMFAIterator, x) = mul!(deepcopy(opt.ZI), opt, x)

# 实现矩阵的伴随算子矩阵向量相乘 
Base.adjoint(opt::MLMFAIterator) = Adjoint(opt)
Base.size(adjopt::Adjoint{T, MLMFAIterator{T, V}}) where {T, V} = size(adjoint(adjopt.parent.Znear))
Base.size(adjopt::Adjoint{T, MLMFAIterator{T, V}}, ind...) where {T, V}  = size(adjoint(adjopt.parent.Znear), ind...)
Base.:*(adjopt::Adjoint{T, MLMFAIterator{T, V}}, x::AbstractVector{S}) where {T, V, S} = mul!(deepcopy(adjopt.parent.ZI), adjopt, x)

include("AdjointIterateOnOctree.jl")

"""
    LinearAlgebra.mul!(y, Zopt::MLMFAIterator, x)

    重载以实现矩阵向量乘积计算
TBW
"""
function LinearAlgebra.mul!(y, adjZopt::Adjoint{T, MLMFAIterator{T, V}}, x) where {T, V}
    Zopt = adjZopt.parent
    # near
    fill!(y, 0)
    mul!(y, adjoint(Zopt.Znear), x)
    # far
    caladjZfarI!(adjZopt, x)
    # cumsum results
    axpy!(1, Zopt.ZI, y)

    return y

end

"""
    LinearAlgebra.mul!(y, Zopt::MLMFAIterator, x)

    重载以实现矩阵向量乘积计算
TBW
"""
function LinearAlgebra.mul!(y::AbstractVector, adjZopt::Adjoint{T, MLMFAIterator{T, V}}, x::AbstractVector, α::Number, β::Number) where {T,V}
    Zopt = adjZopt.parent
    # near adjoint
    mul!(y, adjoint(Zopt.Znear), x, α, β)
    # far
    caladjZfarI!(adjZopt, x)
    # cumsum results
    axpy!(α, Zopt.ZI, y)

    return y

end