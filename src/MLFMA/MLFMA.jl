## 采用多层快速多极子算法加速计算过程
# 载入八叉树信息
include("OctreeInfo.jl")
# 盒子内几何体 id 的计算函数
include("SetGeoIDsInLevelCubes.jl")
# 基于八叉树的迭代器
include("MLMFAIterators.jl")

"""
计算基函数中心的数组，用于方便混合基函数使用时的情况
"""
function getBfsCenter(bfsInfo::Vector{BFT}) where {BFT<:BasisFunctionType}
    
    # 基函数数量
    nbf = length(bfsInfo)
    # 预分配内存基函数公共边中点内存
    bfcenters   =   zeros(eltype(bfsInfo[1].center), (3, nbf))
    # 循环计算
    @threads for ibf in 1:nbf
        bfcenters[:, ibf]   .=  bfsInfo[ibf].center
    end #nbf

    bfcenters
end

"""
计算基函数中心的数组，用于方便混合基函数使用时的情况
"""
function getBfsCenter(bfsInfo::Vector{VT}) where {VT<:AbstractVector}
    # 基函数数量
    nbf = sum(bfs -> length(bfs), bfsInfo)
    # 预分配内存基函数公共边中点内存
    bfcenters   =   zeros(eltype(bfsInfo[1][1].center), (3, nbf))
    # 循环计算
    for bfs in bfsInfo
        @threads for ibf in eachindex(bfs)
            bfcenters[:, bfs[ibf].bfID]   .=  bfs[ibf].center
        end #nbf
    end

    bfcenters
end

"""
三角形面网格直接设置为的叶层盒子边长
"""
function getLeafCubeL(geosInfo::Vector{T}) where {T<:TriangleInfo}
    FT = Precision.FT
    re::FT = MLFMAParams.LEAFCUBESIZE
    re
end

"""
四面体从网格平均尺寸设置整体的叶层盒子边长
"""
function getLeafCubeL(geosInfo::Vector{T}) where {T<:TetrahedraInfo}
    FT = Precision.FT
    re::FT = 7mean(x -> mean(f -> mean(f.edgel), x.faces), geosInfo)/4
    re
end

"""
六面体从网格平均尺寸设置整体的叶层盒子边长
"""
function getLeafCubeL(geosInfo::AbstractVector{T}) where {T<:HexahedraInfo}
    FT = Precision.FT
    re::FT = 2mean(x -> mean(f -> mean(f.edgel), x.faces), geosInfo)
    re
end

"""
混合网格从第2类里的网格平均尺寸设置整体的叶层盒子边长
"""
function getLeafCubeL(geosInfo::Vector{T}) where {T<:AbstractVector}
    getLeafCubeL(geosInfo[end])
end


"""
根据基函数中心位置构建八叉树，并重排基函数信息、将新基函数 ID 赋值给几何元信息数组
返回值：
nLevels::   层数
octree::    得到的八叉树
leafCubeEdgel:: 控制叶层盒子大小
isDistribute:: 控制是否为分布式计算
"""
function getOctreeAndReOrderBFs!(geosInfo, bfsInfo; leafCubeEdgel = getLeafCubeL(geosInfo), nInterp = 6)
    @clock "八叉树构建" begin
        # 所有基函数的中心坐标，用于分组计算
        bfcenters   =  getBfsCenter(bfsInfo)
        
        # 构建八叉树，重新排列基函数
        octree, reOrderID   =   OctreeInfo{eltype(bfcenters), LevelInfo}(bfcenters, leafCubeEdgel; nInterp = nInterp)

        # 叶层ID
        nLevels     =   max(keys(octree.levels)...)
        
        @info "基函数按八叉树重新排序中..."
        # 根据按八叉树重新排序的基函数id重排基函数信息
        reOrderBasisFunctionAndGeoInfo!(reOrderID, geosInfo, bfsInfo)
        # 在叶层找出每个盒子包含的体元id
        setGeoIDsInLeafCubes!(octree.levels[nLevels], bfsInfo)
    end

    return nLevels, octree

end
