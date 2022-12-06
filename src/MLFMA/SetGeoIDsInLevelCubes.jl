"""
用于设置给定层的盒子中包含的几何体，采用 RWG、SWG、RBF 基函数时，八叉树分组依据为基函数，
同一个几何体会被分在不同的基函数上会被分入入不同的盒子，因此邻盒子中的几何体 id 大概率出现重复值。
"""
function setGeoIDsInLeafCubes!(level, bfsInfo::Vector{LBF}) where{LBF<:LinearBasisFunction}
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 叶层盒子数量
    nCubes      =   cubesIndices.stop
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        # 盒子
        cube    =   cubes[iCube]
        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的几何体id
        cubeGeoID       =   zeros(Int, 2*length(cubeBFinterval))
        @inbounds for (i, bf) in enumerate(view(bfsInfo, cubeBFinterval))
            cubeGeoID[(2i-1):(2i)]  .=  bf.inGeo
        end
        # 排序并剔除冗余元素
        unique!(sort!(cubeGeoID))
        # 防止边界元基函数的负部，编号为 0 的盒子出现在索引中
        cubeGeoID[1] == 0 && popfirst!(cubeGeoID)

        cube.geoIDs    =   cubeGeoID
    end

    return nothing
end

"""
用于设置给定层的盒子中包含的几何体，采用 RWG、SWG、RBF 基函数时，八叉树分组依据为基函数，
同一个几何体会被分在不同的基函数上会被分入入不同的盒子，因此邻盒子中的几何体 id 大概率出现重复值。
"""
function setGeoIDsInLeafCubes!(level, bfsInfo::Vector{VT}) where{VT<:AbstractVector}
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 叶层盒子数量
    nCubes      =   cubesIndices.stop
    # 拆分基函数信息
    sbfs    =   bfsInfo[1]
    vbfs    =   bfsInfo[2]

    # 面基函数在所有基函数中的id
    sbfIDs  =   [bf.bfID for bf in sbfs]
    # 体基函数在所有基函数中的id
    vbfIDs  =   [bf.bfID for bf in vbfs]

    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        # 盒子
        cube    =   cubes[iCube]
        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的几何体id
        cubeGeoID       =   zeros(Int, 2*length(cubeBFinterval))
        @inbounds for i in 1:length(cubeBFinterval)
            # 在面、体基函数中？
            n       =   cubeBFinterval[i]
            insbf   =   searchsorted(sbfIDs, n)
            invbf   =   searchsorted(vbfIDs, n)
            !isempty(insbf) && begin
                bf = sbfs[insbf.start]
                cubeGeoID[(2i-1):(2i)]  .=  bf.inGeo
            end
            !isempty(invbf) && begin
                bf = vbfs[invbf.start]
                cubeGeoID[(2i-1):(2i)]  .=  bf.inGeo
            end
        end
        # 排序并剔除冗余元素
        unique!(sort!(cubeGeoID))
        # 防止边界元基函数的负部，编号为 0 的盒子出现在索引中
        cubeGeoID[1] == 0 && popfirst!(cubeGeoID)

        cube.geoIDs    =   cubeGeoID
    end

    return nothing
end


"""
用于设置给定层的盒子中包含的几何体，采用常数基函数时，同一个盒子不会出现重复值。
"""
function setGeoIDsInLeafCubes!(level, bfsInfo::Vector{CBF}) where{CBF<:ConstBasisFunction}
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 叶层盒子数量
    nCubes      =   cubesIndices.stop
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        # 盒子
        cube    =   cubes[iCube]

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的几何体id
        cubeGeoID       =   zeros(Int, length(cubeBFinterval))
        @inbounds for (i, bf) in enumerate(view(bfsInfo, cubeBFinterval))
            cubeGeoID[i]  =  bf.inGeo
        end
        # 排序并剔除冗余元素
        unique!(sort!(cubeGeoID))
        # 写入数据
        cube.geoIDs    =   cubeGeoID
    end

    return nothing
end