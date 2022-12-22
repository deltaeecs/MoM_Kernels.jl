spmul!  =   (@isdefined MKLSparse) ? MKLSparse.BLAS.mul! : SparseArrays.mul!


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
LinearMapType{T}   =   Union{AbstractMatrix{T}, MLMFAIterator{T, VT}} where {T<:Number, VT}

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
    aggSBF, disaggSBF   =   getAggSBFOnLevel(leafLevel, vsCellsInfo, bfsInfo)

    # 给各层的聚合项、解聚项预分配内存
    memoryAllocationOnLevels!(nLevels, levels)

    # 给矩阵向量乘积预分配内存
    ZI  =   zeros(Complex{FT}, nbf)

    return MLMFAIterator{Complex{FT}, Vector}(octree, ZnearCSC, vsCellsInfo, bfsInfo, aggSBF, disaggSBF, ZI)

end

"""
    LinearAlgebra.mul!(y, Zopt::MLMFAIterator, x)

    重载以实现矩阵向量乘积计算
TBW
"""
function LinearAlgebra.mul!(y, Zopt::MLMFAIterator, x)
    # near
    mul!(y, Zopt.Znear, x)
    # far
    calZfarI!(Zopt, x)
    # cumsum results
    axpy!(1, Zopt.ZI, y)

    return y

end

Base.:*(opt::MLMFAIterator, x) = mul!(copy(opt.ZI), opt, x)
