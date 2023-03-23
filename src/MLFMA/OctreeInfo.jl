## 定义八叉树类型以及构造函数
include("LevelInfo.jl")
## 计算相移因子、转移因子的函数
include("PhaseShiftAndTransFactors.jl")

"""
八叉树类
nLevels ::Integer, 叶层ID（定义大盒子为（“0” 层），叶层为第“n”层，nLevels取“n”的值）
leafCubeEdgel::FT，叶层盒子边长
bigCubeLowerCoor::MVec3D{FT}，第0层盒子的角坐标
levels  ::Dict{Int, LevelInfo}，保存各层信息的字典
"""
struct OctreeInfo{FT<:Real, LT<:AbstractLevel}
    nLevels         ::Integer
    leafCubeEdgel   ::FT
    bigCubeLowerCoor::Vec3D{FT}
    levels          ::Dict{IT, LT} where {IT}
end

"""
构建八叉树类
"""
function OctreeInfo{FT, LT}(leafnodes::Matrix{FT}, leafCubeEdgel::FT; nInterp = 6, 
    IPT = get_Interpolation_Method(MLFMAParams.InterpolationMethod)) where{FT<:Real, LT<:AbstractLevel}
    
    println("构造八叉树中...")
    
    # 构建将叶点完全覆盖的“大盒子”, 其中leafCubeEdgelUsed为计算出的实际使用的叶层盒子边长
    nLevels, bigCubeLowerCoor,  leafCubeEdgelUsed  =   setBigCube(leafnodes, leafCubeEdgel)
    nLevels <= 2 && error("层数过少，检查参数！")

    # 创建叶层
    leafLevel       =   LT{Int, FT, IPT{Int, FT}}()
    leafsIDSorted   =   setLevelInfo!(leafLevel, nLevels, leafnodes, leafCubeEdgelUsed, bigCubeLowerCoor)

    # 将叶层写入字典
    levels  =   Dict{Int, LT}(nLevels => leafLevel)
    # 记录按按所在父层盒子排序后的子盒子id
    levelsCubeIDSorted    =   Dict{Int, Vector{Int}}()
    levelsCubeIDSorted[nLevels+1]   =   leafsIDSorted

    # 创建非叶层
    for ilevel in (nLevels - 1):-1:1
        # ilevel的盒子边长
        ilevelCubeEdgel =   leafCubeEdgelUsed*(2^(nLevels - ilevel))
        # 创建父层，同时输出子层盒子排序后的新id
        levels[ilevel]  =   LT{Int, FT, IPT{Int, FT}}()
        levelsCubeIDSorted[ilevel + 1]  =   setLevelInfo!(levels[ilevel], ilevel, levels[ilevel + 1], ilevelCubeEdgel, bigCubeLowerCoor)
    end

    # 从顶层向叶层， 根据排序后的新id重新排列子层盒子，以将同一个父盒子层的盒子相邻排列，这样有利于计算
    reOrderCubeID!(nLevels, levels, levelsCubeIDSorted)
    # 从叶层向顶层，计算每一层的每个盒子包含的基函数的区间
    setBFInterval!(nLevels, levels)
    # 从 (iLevel - 1) 到顶层，计算每层非空盒子包含的非空子盒子在本盒子中的相对位置
    setLevelsCubesKidsIn8!(nLevels, levels)

    # 可以方便地计算远亲组了
    for ilevel in 1:(nLevels - 1)
        setKidLevelFarNeighbors!(levels[ilevel], levels[ilevel + 1])
    end

    println("预计算采样点、插值矩阵、相移因子、转移因子等信息中...")
    # 计算插值信息
    for ilevel in nLevels:-1:2
        levels[ilevel].interpWθϕ    =   interpolationCSCMatCal(levels[ilevel-1].poles, levels[ilevel].poles, nInterp)
    end

    # 预计算层间盒子的相移因子
    setLevelsShiftFactor!(nLevels, levels)

    # 预计算层内转移因子
    setLevelTransFactor!(nLevels, levels)

    println("八叉树构造完毕")
    # 返回八叉树和重排过的叶子（基函数）id
    return OctreeInfo{FT, LT}(nLevels, leafCubeEdgel, bigCubeLowerCoor, levels), leafsIDSorted

end

OctreeInfo(leafnodes::Matrix{FT}, leafCubeEdgel::FT) where{FT<:Real} = OctreeInfo{FT, }(leafnodes::Matrix{FT}, leafCubeEdgel::FT)

"""
实现包含分布式层的
"""
function Base.convert(::Type{OctreeInfo{FT, LTt}}, octree::OctreeInfo{FT, LTo}) where {FT<:Real, LTt<:AbstractLevel, LTo<:AbstractLevel}
    # 不需要转换
    LTo === LTt && return octree
    # 需要转换
    levelst =   Dict{Int, LTt}()
    for (k, v) in octree.levels
        levelst[k] = convert(LTt, v)
    end

    OctreeInfo{FT, LTt}(octree.nLevels, octree.leafCubeEdgel, octree.bigCubeLowerCoor, levelst)

end
"""
计算包围目标的大盒子信息
输入：
nodes::Matrix{FT}，大小为 (3, n) 的用于分割成八叉树的空间点，如基函数的中心坐标，或者为构成网格的所有点
leafCubeEdgel::FT，叶层盒子边长，用于计算总层数和大盒子的坐标信息
"""
function setBigCube(nodes::Matrix{FT}, leafCubeEdgel::FT) where{FT<:Real}

    ## 计算包住目标的盒子的坐标
    # 长、宽、高，附加小距离以完全圈主
    xyzmin  =   minimum(nodes, dims=2)
    xyzmax  =   maximum(nodes, dims=2)
    Δxyz    =   xyzmax .- xyzmin .+ (sqrt(2.0)-1.0)*leafCubeEdgel
    # 包住目标的盒子的边长
    CubeEdgel    =   maximum(Δxyz)
    # 包住目标的盒子的坐标
    CubeCenter   =   SVec3D{FT}(xyzmin .+ xyzmax) ./ 2
    # 总叶层ID
    nLevels         =   ceil(Int, log2(CubeEdgel/leafCubeEdgel))
    
    #  严格使用设定的叶层盒子边长
    ## 根据算出的层数计算第0层的大盒子
    # 第0层盒子的边长
    bigCubeEdgel    =   leafCubeEdgel*2^nLevels
    # 大盒子中心的角坐标
    bigCubeLowerCoor=   CubeCenter .- bigCubeEdgel/2

    return nLevels, bigCubeLowerCoor, leafCubeEdgel
end


"""
根据排序后的新id重新排列子层盒子以及盒子的邻盒子信息，以将同一个父盒子层的盒子相邻排列，这样有利于计算
更新的量：父层盒子的kidsInterval， 本层的盒子顺序，本层盒子的邻盒子id
"""
function reOrderCubeID!(nLevels::Integer, levels::Dict{Int, LV}, levelsCubeIDSorted::Dict{Int, Vector{Int}})  where {LV<:AbstractLevel}
    
    # 从第“2”层开始计算到叶层
    for ilevel in 2:nLevels
        # 先按父层盒子包含子盒子的顺序重新排列cubeIDSorted
        # 父层
        parentLevel =   levels[ilevel - 1]
        # 本层
        level       =   levels[ilevel]
        # 本层原本的 排序后id
        cubeIDSorted=   levelsCubeIDSorted[ilevel]
        # 原地计算新排序
        permute!(cubeIDSorted, vcat([cube.kidsInterval for cube in parentLevel.cubes]...))
        # 更新父层的kidsInterval
        let nCumKidCube = 0
            for pCube in parentLevel.cubes
                # 该父盒子包含的子盒子数
                nkid    =   length(pCube.kidsInterval)
                # 重新计算区间
                pCube.kidsInterval  =   (nCumKidCube + 1) :(nCumKidCube + nkid)
                # 累加
                nCumKidCube +=  nkid
            end #for pCUbe
        end# let
        
        # 按新的cubeIDSorted排序本层盒子
        permute!(level.cubes, cubeIDSorted)
        # 盒子数
        nCubes      =   length(cubeIDSorted)
        # 旧id新id对照
        oldNewIDpair=   hcat(cubeIDSorted, 1:nCubes)
        # 对旧id进行排序得到按旧id排序的新id
        oldNewIDpair=   sortslices(oldNewIDpair, dims = 1, alg = ThreadsX.MergeSort)
        # 再更新邻盒子id
        for cube in level.cubes
            cube.neighbors  .= oldNewIDpair[cube.neighbors, 2]
            sort!(cube.neighbors)
        end #cube
    end # ilevel

    # 对叶层盒子的kidsInterval即基函数区间进行更新
    leafLevel   =   levels[nLevels]
    leafsIDSorted   =   levelsCubeIDSorted[nLevels + 1]
    # 原地计算新排序
    permute!(leafsIDSorted, vcat([cube.kidsInterval for cube in leafLevel.cubes]...))
    # 更新父层的kidsInterval
    let nCumKidCube = 0
        for lCube in leafLevel.cubes
            # 该父盒子包含的子盒子数
            nkid    =   length(lCube.kidsInterval)
            # 重新计算区间
            lCube.kidsInterval  =   (nCumKidCube + 1) :(nCumKidCube + nkid)
            # 累加
            nCumKidCube +=  nkid
        end #for lCUbe
    end# let
    
    return nothing
end #function

"""
根据按八叉树重新排序的id重排基函数信息
"""
function reOrderBasisFunctionAndGeoInfo!(bfReOrderID::Vector{IT}, geosInfo::Vector{VCellT}, bfsInfo::Vector{BFT}) where {IT<:Integer, BFT<:BasisFunctionType, VCellT<:VSCellType}
    
    # 重新计算旧基函数id到新的基函数id的映射
    # 旧id新id对照
    oldNewIDpair    =   hcat(bfReOrderID, 1:length(bfReOrderID))
    # 对旧id进行排序得到按旧id排序的新id
    oldNewIDpair    =   sortslices(oldNewIDpair, dims = 1, alg = ThreadsX.MergeSort)
    # 更新基函数ID 并按其排序
    @threads for bf in bfsInfo
        bf.bfID = oldNewIDpair[bf.bfID, 2]
    end
    sort!(bfsInfo, by = x -> x.bfID, alg = ThreadsX.MergeSort)

    # 再更新几何单元内的基函数 id
    @threads for geo in geosInfo
        for ii in eachindex(geo.inBfsID)
            geo.inBfsID[ii]  = oldNewIDpair[geo.inBfsID[ii], 2]
        end
    end # cube

    return nothing
end

"""
根据按八叉树重新排序的id重排基函数信息
"""
function reOrderBasisFunctionAndGeoInfo!(bfReOrderID::Vector{IT}, geosInfo::Vector{VT1}, bfsInfo::Vector{VT2}) where {IT<:Integer, VT1<:AbstractVector, VT2<:AbstractVector}

    # 其次重新计算旧基函数id到新的基函数id的映射
    # 旧id新id对照
    oldNewIDpair    =   hcat(bfReOrderID, 1:length(bfReOrderID))
    # 对旧id进行排序得到按旧id排序的新id
    oldNewIDpair    =   sortslices(oldNewIDpair, dims = 1, alg = ThreadsX.MergeSort)

    for bfs in bfsInfo
        # 更新基函数ID 并按其排序
        for bf in bfs
            bf.bfID = oldNewIDpair[bf.bfID, 2]
        end
        sort!(bfs, by = x -> x.bfID, alg = ThreadsX.MergeSort)
    end
    # 再更新几何单元内的基函数 id
    for geos in geosInfo
        istri = eltype(geos)<:TriangleInfo
        @threads for vsCellInfo in geos
            for ii in eachindex(vsCellInfo.inBfsID)
                istri && vsCellInfo.inBfsID[ii] == 0 && continue
                vsCellInfo.inBfsID[ii]  = oldNewIDpair[vsCellInfo.inBfsID[ii], 2]
            end
        end # cube
    end

    return nothing
end


"""
根据按八叉树重新排序的id重排基函数信息，此函数适用于 PWC 基函数的情况
"""
function reOrderBasisFunctionAndGeoInfo!(reOrderID::Vector{IT}, vsCellsInfo::Vector{VCellT}, ::Val{:PWC}) where {IT<:Integer, VCellT<:VSCellType}

    # 首先重排rwg基函数信息
    permute!(vsCellsInfo, reOrderID)
    # 再更新几何单元内的基函数 id
    for vsi in 1:length(vsCellsInfo)
        vsCellsInfo[vsi].inBfsID    .=  (3(vsi-1)+1):(3vsi)
    end #cube

    return nothing
end

"""
根据已经排序好的层的盒子信息，从叶层到顶层更新盒子包含的基函数区间
"""
function setBFInterval!(nLevels::Integer, levels::Dict{Int, LV}) where {LV<:AbstractLevel}

    # 先更新子层盒子，与子层盒子的 kidsInterval属性相同
    leafCubes=  levels[nLevels].cubes
    for cube in leafCubes
        cube.bfInterval = cube.kidsInterval
    end
    
    # 从第“2”层开始计算到叶层
    for ilevel in nLevels-1:-1:1
        # 本层盒子
        tlevelCubes =   levels[ilevel].cubes
        # 子层盒子
        klevelCubes =   levels[ilevel+1].cubes
        # 子层盒子包含的基函数
        for tlevelCube in tlevelCubes
            # 本盒子的子盒子
            tklevelCubes    =   view(klevelCubes, tlevelCube.kidsInterval)

            # 子盒子的基函数起点
            kidsbfstart =   tklevelCubes[1].bfInterval[1]
            # 子盒子的基函数终点
            kidsbfend   =   tklevelCubes[end].bfInterval[end]
            # 写入本盒子的基函数区间
            tlevelCube.bfInterval   =   kidsbfstart:kidsbfend
        end #tlevelCube
    end #ilevel

end #function

"""
计算（nLevel-1）-2 层每层的非空盒子的非空子盒子在其8个子盒子中的位置
"""
function setLevelsCubesKidsIn8!(nLevels::Integer, levels::Dict{Int, LV}) where {LV<:AbstractLevel}
    for iLevel in (nLevels - 1):-1:2
        # 本层盒子
        tCubes  =   levels[iLevel].cubes
        # 子层盒子
        kCubes  =   levels[iLevel + 1].cubes
        for iCube in eachindex(tCubes)
            # 本盒子
            cube    =   tCubes[iCube]
            # 本盒子的3Did * 2
            cubeID3D4   =   cube.ID3D .* 4 .- 2
            # 对子盒子循环
            for (ii, kCube) in enumerate(view(kCubes, cube.kidsInterval))
                # 子盒子的 3Did 与 cubeID3D4 的差值为子盒子相对本盒子的位置
                relativeOffset  =   (kCube.ID3D .* 2 .- 1) .- cubeID3D4
                relativeOffsetTrunc =   [relativeOffset[i] == -1 ? 0 : 1 for i in 1:3]
                relativeOffset1D    =   1 + relativeOffsetTrunc[1] + relativeOffsetTrunc[2]*2 + relativeOffsetTrunc[3]*4
                # 写入数据
                cube.kidsIn8[ii]    =   relativeOffset1D
            end # (ii, kCube)
        end # iCube
    end #iLevel
    return nothing
end

##
##