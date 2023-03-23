## 本文件包含八叉树的层类文件

include("IntegralInterpolationInfo.jl")

"""
盒子信息，包括
子层盒子id的区间、
包含的基函数区间、
非空子盒子在8个子盒子中的id等、
包含的网格如三角形、四面体的id，以基函数进行分，因此边界上的同一个网格可能被分到不同的盒子内。
邻盒子的id、远亲盒子的id、本盒子在本层的三维整数坐标、本盒子在本层的三维全局坐标
"""
mutable struct CubeInfo{IT<:Integer, FT<:Real}
    kidsInterval    ::UnitRange{IT}
    bfInterval      ::UnitRange{IT}
    kidsIn8         ::Vector{IT}
    geoIDs          ::Vector{IT}
    neighbors       ::Vector{IT}
    farneighbors    ::Vector{IT}
    ID3D            ::MVec3D{IT}
    center          ::MVec3D{FT}
end


abstract type AbstractLevel <: Any end

"""
层信息
ID          ::IT，层序号
L           ::IT, 本层截断项数
cubes       ::Vector{CubeInfo{IT, FT}} 包含每一个盒子信息的向量
cubeEdgel   ::FT，本层盒子的边长
poles       ::PolesInfo{IT, FT}, 多极子采样信息
interpWθϕ   ::InterpInfo{IT, FT}, 插值信息
aggS        ::Array{Complex{FT}, 3}， 聚合项
disaggG     ::Array{Complex{FT}, 3}， 解聚项
phaseShift2Kids  ::Array{Complex{FT}, 3}，本层盒子到子层盒子的相移因子 
αTrans      ::Array{Complex{FT}, 3}， 本层盒子远亲组之间的转移因子，根据相对位置共有 7^3 - 3^3 = 316 个
αTransIndex ::Array{IT, 2}, 远亲盒子的相对位置到其转移因子在所有转移因子数组的索引
"""
mutable struct LevelInfo{IT<:Integer, FT<:Real, IPT} <: AbstractLevel
    ID          ::IT
    L           ::IT
    nCubes      ::IT
    cubes       ::Vector{CubeInfo{IT, FT}}
    cubeEdgel   ::FT
    poles       ::PolesInfo{FT}
    interpWθϕ   ::IPT
    aggS        ::Array{Complex{FT}, 3}
    disaggG     ::Array{Complex{FT}, 3}
    phaseShift2Kids     ::Array{Complex{FT}, 2}
    phaseShiftFromKids  ::Array{Complex{FT}, 2}
    αTrans      ::Array{Complex{FT}, 2}
    αTransIndex ::OffsetArray{IT, 3, Array{IT, 3}}

    LevelInfo{IT, FT, IPT}() where {IT<:Integer, FT<:Real, IPT<:InterpInfo} = new{IT, FT, IPT}()
    LevelInfo{IT, FT, IPT}(  ID, L, nCubes, cubes, cubeEdgel, poles, interpWθϕ, aggS, disaggG,
                        phaseShift2Kids, phaseShiftFromKids, αTrans,  αTransIndex) where {IT<:Integer, FT<:Real, IPT<:InterpInfo} = 
            new{IT, FT, IPT}(ID, L, nCubes, cubes, cubeEdgel, poles, interpWθϕ, aggS, disaggG,
                        phaseShift2Kids, phaseShiftFromKids, αTrans,  αTransIndex)
end


"""
叶层LevelInfo的构造函数，输入为空间三维坐标数组
nLevels::IT，层 数，亦为叶层层ID
leafnodes::Matrix{FT},大小为 (3, n) 的用于分割成八叉树的空间点，如基函数的中心坐标
cubeEdgel::FT，叶层盒子边长
bigCubeLowerCoor::Vec3D{FT}， 大盒子的角坐标
"""
function setLevelInfo!(level::LT, nLevels::Integer, leafnodes::Matrix{FT},
    cubeEdgel::FT, bigCubeLowerCoor::Vec3D{FT}) where{LT<:AbstractLevel, FT<:Real}
    # 计算
    nleaves =   size(leafnodes, 2)
    # 每一个节点所在盒子的3Did（3列）+节点编号（1列）
    nodesInCubeID3D =   zeros(Int, (nleaves, 4))
    nodesInCubeID3D[:, 4]   =   1:nleaves
    @threads for ileaf in 1:nleaves
        nodesInCubeID3D[ileaf, 1:3]    .=  ceil.(Int, (leafnodes[:, ileaf] .- bigCubeLowerCoor) ./ cubeEdgel)
    end
    
    # 通过排序，将相同盒子的基函数放在一起
    nodesInCubeID3D    .=   sortslices(nodesInCubeID3D, dims = 1, alg = ThreadsX.MergeSort)
    # 排序后的基函数信息（用于重新排列基函数编号，将相同盒子的基函数放在一起，便于计算）
    kidsSorted  =   nodesInCubeID3D[:, 4]

    # 计算盒子包含的基函数信息，先考虑最坏情况分配足够大的数组
    kidsIntervals   =   ones(Int, nleaves)
    let icube = 1
        for ileaf in 1:(nleaves-1)
            # 判断与后一个是否为同一个
            if nodesInCubeID3D[ileaf, 1:3] != nodesInCubeID3D[ileaf+1, 1:3]
                icube   +=    1
                kidsIntervals[icube]    =   ileaf+1
            end
        end
        # 把基函数数量+1放到数组最后方便计算
        icube  +=    1
        kidsIntervals[icube]    =   nleaves+1
        # 将kidsIntervals恢复到(盒子数+1)大
        resize!(kidsIntervals, icube)
    end#let
    # 本层非空盒子数目
    nCubes  =   length(kidsIntervals) - 1
    # 创建盒子包含的基函数区间切片向量并计算
    kidsSlice       =   [kidsIntervals[i]:(kidsIntervals[i+1]-1) for i in 1:nCubes]
    # 盒子的三维id
    cubesID3D       =   nodesInCubeID3D[kidsIntervals[1:(end-1)], 1:3]

    ## 开始查找每一个盒子的邻盒子id
    cubesNeighbors  =   searchNearCubes(cubesID3D, nLevels)
    
    # 构建数组保存盒子信息
    cubesInfo   =   Vector{CubeInfo{Int, FT}}(undef, nCubes)
    @threads for icube in 1:nCubes
        # 盒子中心
        center  =   bigCubeLowerCoor .+ (cubesID3D[icube, :] .- 1/2) .* cubeEdgel
        cubesInfo[icube] =   CubeInfo{Int, FT}(kidsSlice[icube],  0:0, Int[], Int[], cubesNeighbors[icube], Int[], cubesID3D[icube, :], center)
    end #for
    
    # 计算截断项数和角谱空间采样多极子信息
    L::Int, poles    =   levelIntegralInfoCal(cubeEdgel, Val(MLFMAParams.InterpolationMethod))

    # 将相关项写入level
    level.ID        =   nLevels
    level.L         =   L
    level.nCubes    =   nCubes
    level.cubes     =   cubesInfo
    level.cubeEdgel =   cubeEdgel
    level.poles     =   poles


    # 层按盒子排序后的id
    return kidsSorted

end

"""
非叶层LevelInfo的构造函数，输入为空间三维坐标数组
levelID::计算层的id
leafnodes::Matrix{FT},大小为 (3, n) 的用于分割成八叉树的空间点，如基函数的中心坐标
cubeEdgel::FT，本层盒子边长
"""
function setLevelInfo!(level, levelID::Integer, kidLevel, cubeEdgel::FT, bigCubeLowerCoor::Vec3D{FT}) where{FT<:Real}
    # 计算
    nkidCubes   =   length(kidLevel.cubes)
    # 每一个子盒子所在父盒子的3Did（3列）+子盒子点编号（1列）
    kidCubesInCubeID3D          =   zeros(Int, (nkidCubes, 4))
    kidCubesInCubeID3D[:, 4]   .=   1:nkidCubes
    @threads for ikidCube in 1:nkidCubes
        kidCubesInCubeID3D[ikidCube, 1:3]    .=  ceil.(Int, kidLevel.cubes[ikidCube].ID3D ./ 2)
    end
    
    # 通过排序，将相同父盒子的子盒子放在一起
    kidCubesInCubeID3D .=   sortslices(kidCubesInCubeID3D, dims = 1, alg = ThreadsX.MergeSort)
    # 排序后的子盒子id信息（用于重新排列子盒子编号，将相同父盒子的子盒子放在一起，便于计算）
    kidCubesSorted      =   kidCubesInCubeID3D[:, 4]

    # 计算盒子包含的子盒子信息，先考虑最坏情况分配足够大的数组
    kidsIntervals   =   ones(Int, nkidCubes+1)
    let icube = 1
        for ikidCube in 1:(nkidCubes - 1)
            # 判断与后一个是否为同一个
            if kidCubesInCubeID3D[ikidCube, 1:3] != kidCubesInCubeID3D[ikidCube+1, 1:3]
                icube   +=    1
                kidsIntervals[icube]    =   ikidCube + 1
            end
        end
        # 把子盒子数量+1放到数组最后方便计算
        icube  +=    1
        kidsIntervals[icube]    =   nkidCubes+1
        # 将kidsIntervals恢复到(盒子数+1)大
        resize!(kidsIntervals, icube)
    end#let
    # 本层非空盒子数目
    nCubes  =   length(kidsIntervals) - 1
    # 创建盒子包含的子盒子区间切片向量并计算
    kidsSlice       =   [kidsIntervals[i]:(kidsIntervals[i+1]-1) for i in 1:nCubes]
    kidsIn8 =   [Vector{Int}(undef, length(kidSlice)) for kidSlice in kidsSlice]
    
    # 盒子的三维id
    cubesID3D       =   kidCubesInCubeID3D[kidsIntervals[1:(end-1)], 1:3]

    ## 开始查找每一个盒子的邻盒子id
    cubesNeighbors  =   searchNearCubes(cubesID3D, levelID)

    # 构建结构化数组保存盒子信息
    cubesInfo   =   Vector{CubeInfo{Int, FT}}(undef, nCubes)
    for icube in 1:nCubes
        # 盒子中心
        center  =   bigCubeLowerCoor .+ (cubesID3D[icube, :] .- 1/2) .* cubeEdgel
        cubesInfo[icube] =   CubeInfo{Int, FT}(kidsSlice[icube], 0:0, kidsIn8[icube], Int[], cubesNeighbors[icube], Int[], cubesID3D[icube, :], center)
    end #for

    # 计算截断项数和角谱空间采样多极子信息
    L::Int, poles    =   levelIntegralInfoCal(cubeEdgel, Val(MLFMAParams.InterpolationMethod))

    # 将相关项写入level
    level.ID        =   levelID
    level.L         =   L
    level.nCubes    =   nCubes
    level.cubes     =   cubesInfo
    level.cubeEdgel =   cubeEdgel
    level.poles     =   poles

    # 返回子层按盒子排序后的id
    return kidCubesSorted

end

"""
用于寻找邻盒子的函数
输入
cubesID3D::Matrix{Int}，(n×3)盒子在本层的三维坐标
levelID::Integer       层编号，（定义大盒子为（“0” 层），叶层为第“n”层
"""
function searchNearCubes(cubesID3D::Matrix{IT}, levelID::Integer) where {IT<:Integer}
    # 盒子数
    nCubes      =   size(cubesID3D, 1)
    
    # 大盒子能容纳该层盒子的个数
    maxCubes1D  =   2^levelID
    cubesID1D   =   IT[((cubeID3D[1] - 1)*maxCubes1D^2 + (cubeID3D[2]-1)*maxCubes1D + cubeID3D[3]) for cubeID3D in eachrow(cubesID3D)]
    
    # 预分配内存
    neighbors   =   Vector{Vector{IT}}(undef, nCubes)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        # 该盒子的3DID
        cubeID3D    =    cubesID3D[iCube, :]

        # 三个维度上邻盒子的偏置（若盒子在中心都为(-1:1)，在边缘则需要去除超出边界部分）
        neighborsOffsets    =   UnitRange{IT}[-1:1, -1:1, -1:1]
        # 判断是否为边界，进行截断
        for ii in 1:3
            cubeID3D[ii] == 1           && (neighborsOffsets[ii] =  0:1; continue)
            cubeID3D[ii] == maxCubes1D  && (neighborsOffsets[ii] = -1:0)
        end # ii
        # 找出所有的邻盒子并判断保存非空的
        neighborsNonEmpty   =   Vector{IT}(undef, 27)
        nonEmptyNum         =   0
        for offsetx in neighborsOffsets[1]
            for offsety in neighborsOffsets[2]
                for offsetz in neighborsOffsets[3]
                    # 是自身吗
                    ((offsetx == 0) & (offsety == 0) & (offsetz == 0)) && begin
                        nonEmptyNum     += 1
                        neighborsNonEmpty[nonEmptyNum]  =   iCube
                        continue
                    end #begin
                    # 对应邻盒子的一维id
                    neighbor1DID    =   (cubeID3D[1] - 1 + offsetx)*maxCubes1D^2 + (cubeID3D[2] - 1 + offsety)*maxCubes1D + cubeID3D[3] + offsetz
                    # 查找对应位置，无匹配则返回一个
                    isNeighborEmpty =   searchsorted(cubesID1D, neighbor1DID)
                    # 这个id非空吗？，是则保存
                    ~isempty(isNeighborEmpty) && begin
                        nonEmptyNum     += 1
                        neighborsNonEmpty[nonEmptyNum]  =   isNeighborEmpty[1]
                    end #begin
                end #for offsetz
            end #for offsety
        end #for offsetx
        # 根据非空邻盒子实际大小截断向量
        resize!(neighborsNonEmpty, nonEmptyNum)
        # 排序
        sort!(neighborsNonEmpty)
        # 保存
        neighbors[iCube] =  neighborsNonEmpty
    end #iCube
    # 返回
    return neighbors
end
##


"""
寻找子层的远亲盒子
输入::
thisLevel::LevelInfo{IT, FT, IPT}, 本层信息
kidLevel::LevelInfo{IT, FT, IPT}， 子层信息
"""
function setKidLevelFarNeighbors!(thisLevel, kidLevel)
    # 所有本层盒子
    cubes   =   thisLevel.cubes
    # 所有子层盒子
    kidCubes=   kidLevel.cubes
    # 对本层盒子循环
    @threads for icube in eachindex(cubes)
        cube        =   cubes[icube]
        # 本盒子的邻盒子
        neighbors   =   cube.neighbors
        # 本盒子的子盒子ID
        kidCubesID  =   cube.kidsInterval

        # 初始化邻盒子的子盒子
        neighborsKids   =   Int[]
        for iNeighbor in neighbors
            iNeighbor !=  icube && push!(neighborsKids, cubes[iNeighbor].kidsInterval...)
        end #iNeighbor

        # 对子盒子循环
        for iKidCube in kidCubesID
            # 子盒子
            kidCube     =   kidCubes[iKidCube]
            # 子盒子的邻盒子
            kidCubeNeighbors    =    kidCube.neighbors
            # 写入子盒子的远亲
            kidCube.farneighbors  =   sort!(setdiff(neighborsKids, kidCubeNeighbors))
        end #iKidCube
    end #iCube

    return nothing
end


"""
预分配各层上的聚合项、解聚项
"""
function memoryAllocationOnLevels!(nLevels::Integer, levels::Dict{IT, LV}) where{IT<:Integer, LV<:LevelInfo}

    # 用到的浮点数类型
    FT = typeof(levels[nLevels].cubeEdgel) 

    for iLevel in nLevels:-1:2
        level   =   levels[iLevel]
        # 本层盒子数
        nCubes  =   level.nCubes
        # 多极子数
        sizePoles   =   length(level.poles.r̂sθsϕs)

        # 开始预分配内存
        # 聚合项
        aggS    =   zeros(Complex{FT}, (sizePoles, 2, nCubes))
        # 解聚项
        disaggG =   zeros(Complex{FT}, (sizePoles, 2, nCubes))
        # 保存
        level.aggS  =   aggS
        level.disaggG = disaggG
    end

    return nothing

end