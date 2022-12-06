"""
采用 RWG 基函数计算 EFIE 面积分（SIE）阻抗矩阵近场元并将结果放在Znear稀疏矩阵中
"""
function calZnearChunkEFIEonCube!(iCube::Int, cubes, 
    geosInfo::AbstractVector{TriangleInfo{IT, FT}},
    Znear, ::Type{BFT}) where {IT<:Integer, FT<:Real, BFT<:LinearBasisFunction}
    
    # 本盒子信息
    cube    =   cubes[iCube]
    # 常数
    Rsglr       =   Params.Rsglr

    # 盒子里的基函数区间
    cubeBFinterval  =   cube.bfInterval
    # 找出对应的三角形id
    cubeGeoID       =   cube.geoIDs

    # 找出所有邻盒子包含的基函数总数、id
    # nNearCubeBFs        =   0 
    # nearCubeBFinterval  =   Vector{UnitRange{IT}}(undef, length(cube.neighbors))
    # 邻盒子中的几何体数
    nNearCubeGeos       =   0 
    # 实际用到的邻盒子数
    nNearCubeCal    =   0
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        nNearCubeCal += 1
        # 邻盒子
        nearCube    =   cubes[jNearCube]

        # # 盒子里的基函数索引
        # nearCubeBFinterval[nNearCubeCal]   =   nearCube.bfInterval
        # # 累加基函数数量
        # nNearCubeBFs    +=  length(nearCube.bfInterval)
        # 累加计算几何体数量
        nNearCubeGeos   +=  length(nearCube.geoIDs)
        
    end
    # nearCubeBFinterval 是按所有邻盒子大小初始化的，此处要调整
    # resize!(nearCubeBFinterval, nNearCubeCal)
    # 收集所有邻盒子（不包括 id 更大的）内的基函数
    # nearCubeBFindices::Vector{IT}   =   vcat(nearCubeBFinterval...)
    # # 排序以便后面更好地使用
    # sort!(nearCubeBFindices)
    nearCubeBFindices = Znear.colIndices

    # 找出对应的三角形id
    nearCubesGeoID  =   zeros(IT, nNearCubeGeos)
    nNearCubeGeoPtr =   1
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        # 邻盒子
        nearCube    =   cubes[jNearCube]
        # 写入邻三角形
        nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs)-1)]  .=   nearCube.geoIDs
        nNearCubeGeoPtr  +=   length(nearCube.geoIDs)
    end
    # 排序并剔除冗余元素
    unique!(sort!(nearCubesGeoID))
    # fetch 到本地
    geosInfolw  =   fetchDVec2LW(geosInfo, nearCubesGeoID)

    # 对盒子内三角形循环
    # @inbounds for iGeo in eachindex(view(geosInfo, cubeGeoID))
    @inbounds for tid in cubeGeoID
        # 局域的场三角形
        # tid   =   cubeGeoID[iGeo]
        trit  =   geosInfolw[tid]
        #= 场三角形与源三角形在不在一个盒子？因为程序利用了PEC目标的CFIE矩阵的对称性
        进行对称位置阻抗矩阵元的计算，要避免对同一个盒子内阻抗矩阵元的重复计算 =#
        # tins  =   nearCubesGeoID[iGeo] in cubeGeoID
        # tins  =   !isempty(searchsorted(cubeGeoID, nearCubesGeoID[iGeo]))
        # 测试三角形包含的三个测试基函数是否在当前盒子（测试盒子）的基函数（测试基函数）区间内
        msInInterval    =   [m in cubeBFinterval for m in trit.inBfsID]

        for sid in nearCubesGeoID
            # 源三角形
            tris    =   geosInfolw[sid]
            # 场源距离
            Rts     =   dist(trit.center, tris.center)
            # 源三角形包含的三个源基函数是否在所有邻盒子（源盒子）的基函数（源基函数）区间内
            nsInInterval    =   [!isempty(searchsorted(nearCubeBFindices, n)) for n in tris.inBfsID]
            # 判断二者远近，调用不同精度的矩阵元处理函数
            if tid == sid
                # 计算三角形相关的(3*3)个矩阵元的结果
                Zts  =  EFIEOnTris(trit)
                # 写入数据
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = tris.inBfsID[ni]
                    # 判断边是不是基函数（边缘不构建半基函数时适用）
                    (m == 0 || n == 0) && continue
                    # 往矩阵填充结果
                    # 判断是不是在源盒子、场盒子包含的区间内
                    ((msInInterval[mi] && nsInInterval[ni])) && begin
                        Znear[m, n] += Zts[mi, ni]
                    end
                end
            elseif Rts < Rsglr
                # 需要进行近奇异性处理的场源三角形
                Zts    =   EFIEOnNearTris(trit, tris)
                # 写入数据
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = tris.inBfsID[ni]

                    # 判断边是不是基函数（边缘不算）
                    (m == 0 || n == 0) && continue
                    ## 分布式避免数据通信不再利用对称性填充
                    # 避免线程锁的矩阵元循环方式下产生的条件
                    # (tid > sid) && (m in cubeBFinterval) && continue
                    # 判断是不是在源盒子、场盒子包含的区间内
                    ((msInInterval[mi] & nsInInterval[ni])) && begin
                        Znear[m, n] += Zts[mi, ni]
                    end
                end
            else
                # 正常高斯求积
                # 计算三角形相关的(3*3)个矩阵元的结果
                Zts    =   EFIEOnTris(trit, tris)
                
                # 写入数据
                
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = tris.inBfsID[ni]

                    # 判断边是不是基函数（边缘不算）
                    (m == 0 || n == 0) && continue
                    # 避免线程锁的矩阵元循环方式下产生的条件
                    (tid > sid) && (m in cubeBFinterval) && continue

                    # 判断是不是在源盒子、场盒子包含的区间内
                    ((msInInterval[mi] & nsInInterval[ni])) && begin
                        Znear[m, n] += Zts[mi, ni]
                    end
                end
                
            end # if
        end #jGeo
    end #iGeo

    return nothing
end

"""
采用 SWG 基函数计算网格元 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在Znear稀疏矩阵中
"""
function calZnearChunkEFIEonCube!(iCube::Int, cubes, 
    geosInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}},
    Znear, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:LinearBasisFunction}

    # 本盒子信息
    cube    =   cubes[iCube]
    # 常数
    Rsglr       =   Params.Rsglr

    # 盒子里的基函数区间
    cubeBFinterval  =   cube.bfInterval
    # 找出对应的三角形id
    cubeGeoID       =   cube.geoIDs

    # 找出所有邻盒子包含的基函数总数、id
    # nNearCubeBFs        =   0 
    # nearCubeBFinterval  =   Vector{UnitRange{IT}}(undef, length(cube.neighbors))
    # 邻盒子中的几何体数
    nNearCubeGeos       =   0 
    # 实际用到的邻盒子数
    nNearCubeCal    =   0
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        nNearCubeCal += 1
        # 邻盒子
        nearCube    =   cubes[jNearCube]

        # # 盒子里的基函数索引
        # nearCubeBFinterval[nNearCubeCal]   =   nearCube.bfInterval
        # # 累加基函数数量
        # nNearCubeBFs    +=  length(nearCube.bfInterval)
        # 累加计算几何体数量
        nNearCubeGeos   +=  length(nearCube.geoIDs)
        
    end
    # nearCubeBFinterval 是按所有邻盒子大小初始化的，此处要调整
    # resize!(nearCubeBFinterval, nNearCubeCal)
    # 收集所有邻盒子（不包括 id 更大的）内的基函数
    # nearCubeBFindices::Vector{IT}   =   vcat(nearCubeBFinterval...)
    # # 排序以便后面更好地使用
    # sort!(nearCubeBFindices)
    nearCubeBFindices = Znear.colIndices

    # 找出对应的三角形id
    nearCubesGeoID  =   zeros(IT, nNearCubeGeos)
    nNearCubeGeoPtr =   1
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        # 邻盒子
        nearCube    =   cubes[jNearCube]
        # 写入邻三角形
        nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs)-1)]  .=   nearCube.geoIDs
        nNearCubeGeoPtr  +=   length(nearCube.geoIDs)
    end
    # 排序并剔除冗余元素
    unique!(sort!(nearCubesGeoID))
    # fetch 到本地
    geosInfolw  =   fetchDVec2LW(geosInfo, nearCubesGeoID)
    # 是否为偏置数组（用于混合网格）
    geoInterval =   getGeosInterval(geosInfo)

    # 对场盒子内四面体循环
    # @inbounds for iGeo in 1:length(nearCubesGeoID)
    @inbounds for tid in cubeGeoID
        
        # 局域的场四面体
        !(tid in geoInterval) && continue
        geot  =   geosInfolw[tid]
        #= 场四面体与源四面体在不在一个盒子？因为程序利用了目标的EFIE矩阵的对称性
        进行对称位置阻抗矩阵元的计算，要避免对同一个盒子内阻抗矩阵元的重复计算 =#
        # tins  =   nearCubesGeoID[iGeo] in cubeGeoID
        # tins  =   !isempty(searchsorted(cubeGeoID, nearCubesGeoID[iGeo]))
        # 测试四面体包含的四个测试基函数是否在当前盒子（测试盒子）的基函数（测试基函数）区间内
        msInInterval    =   [m in cubeBFinterval for m in geot.inBfsID]
        # 局部判断奇异性距离
        Rsglrlc =   Rsglr/sqrt(norm(geot.ε)/ε_0)
        # 对源四面体循环
        for sid in nearCubesGeoID
            
            !(sid in geoInterval) && continue
            # 源四面体
            geos    =   geosInfolw[sid]
            # 场源距离
            Rts     =   dist(geot.center, geos.center)
            # 源四面体包含的四个源基函数是否在所有邻盒子（源盒子）的基函数（源基函数）区间内
            nsInInterval    =   [!isempty(searchsorted(nearCubeBFindices, n)) for n in geos.inBfsID]

            # 判断二者远近，调用不同精度的矩阵元处理函数
            if tid == sid
                # 重合场源四面体
                Zts     =   EFIEOnTetraSWG(geot)
                # 写入数据
                for ni in 1:4, mi in 1:4
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geot.inBfsID[ni]
                    # 往矩阵填充结果
                    # 判断是不是在源盒子、场盒子包含的区间内
                    (msInInterval[mi] && nsInInterval[ni]) && begin
                        Znear[m, n] += Zts[mi, ni]
                    end # begin
                end

            elseif Rts < Rsglrlc
                # 需要进行近奇异性处理的场源四面体
                Zts, _    =   EFIEOnNearTetrasSWG(geot, geos)
                # 写入数据
                for ni in 1:4, mi in 1:4
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geos.inBfsID[ni]

                    # 避免线程锁的矩阵元循环方式下产生的条件
                    (tid > sid) && (m in cubeBFinterval) && continue

                    # 判断是不是在源盒子、场盒子包含的区间内
                    ((msInInterval[mi] && nsInInterval[ni])) && begin
                        Znear[m, n] += Zts[mi, ni]
                    end
                end
            else
                # 正常高斯求积
                # 计算四面体相关的(4*4)个矩阵元的结果
                Zts, _    =   EFIEOnTetrasSWG(geot, geos)
                # 写入数据
                for ni in 1:4, mi in 1:4
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geos.inBfsID[ni]
                    # 避免线程锁的矩阵元循环方式下产生的条件
                    (tid > sid) && (m in cubeBFinterval) && continue

                    # 判断是不是在源盒子、场盒子包含的区间内
                    ((msInInterval[mi] && nsInInterval[ni])) && begin
                        Znear[m, n] += Zts[mi, ni]
                    end
                end
                
            end # if
        end #jGeo
    end #iGeo

    return nothing
end



"""
采用 PWC 基函数计算四面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在Znear稀疏矩阵中
"""
function calZnearChunkEFIEonCube!(iCube::Int, cubes, 
    geosInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}},
    Znear, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:ConstBasisFunction}
    
    # 本盒子信息
    cube    =   cubes[iCube]
    # 常数
    Rsglr       =   Params.Rsglr

    # 找出对应的三角形id
    cubeGeoID       =   cube.geoIDs

    discreteJ::Bool =   (SimulationParams.discreteVar === "J")

    # 找出所有邻盒子包含的基函数总数、id
    # nNearCubeBFs        =   0 
    # nearCubeBFinterval  =   Vector{UnitRange{IT}}(undef, length(cube.neighbors))
    # 邻盒子中的几何体数
    nNearCubeGeos       =   0 
    # 实际用到的邻盒子数
    nNearCubeCal    =   0
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        nNearCubeCal += 1
        # 邻盒子
        nearCube    =   cubes[jNearCube]

        # # 盒子里的基函数索引
        # nearCubeBFinterval[nNearCubeCal]   =   nearCube.bfInterval
        # # 累加基函数数量
        # nNearCubeBFs    +=  length(nearCube.bfInterval)
        # 累加计算几何体数量
        nNearCubeGeos   +=  length(nearCube.geoIDs)
        
    end
    # nearCubeBFinterval 是按所有邻盒子大小初始化的，此处要调整
    # resize!(nearCubeBFinterval, nNearCubeCal)
    # 收集所有邻盒子（不包括 id 更大的）内的基函数
    # nearCubeBFindices::Vector{IT}   =   vcat(nearCubeBFinterval...)
    # # 排序以便后面更好地使用
    # sort!(nearCubeBFindices)
    nearCubeBFindices = Znear.colIndices

    # 找出对应的三角形id
    nearCubesGeoID  =   zeros(IT, nNearCubeGeos)
    nNearCubeGeoPtr =   1
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        # 邻盒子
        nearCube    =   cubes[jNearCube]
        # 写入邻三角形
        nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs)-1)]  .=   nearCube.geoIDs
        nNearCubeGeoPtr  +=   length(nearCube.geoIDs)
    end
    # 排序并剔除冗余元素
    unique!(sort!(nearCubesGeoID))
    # fetch 到本地
    geosInfolw  =   fetchDVec2LW(geosInfo, nearCubesGeoID)
    # 是否为偏置数组（用于混合网格）
    geoInterval =   getGeosInterval(geosInfo)

    # 对场盒子内四面体循环
    # @inbounds for iGeo in 1:length(nearCubesGeoID)
    @inbounds for tid in cubeGeoID
        
        # 局域的场四面体
        !(tid in geoInterval) && continue
        geot  =   geosInfolw[tid]
        # 场四面体介质对比度
        κₜ  =   geot.κ
        #= 场四面体与源四面体在不在一个盒子？因为程序利用了目标的EFIE矩阵的对称性
        进行对称位置阻抗矩阵元的计算，要避免对同一个盒子内阻抗矩阵元的重复计算 =#
        # tins  =   nearCubesGeoID[iGeo] in cubeGeoID
        # tins  =   !isempty(searchsorted(cubeGeoID, nearCubesGeoID[iGeo]))
        # 局部判断奇异性距离
        Rsglrlc =   Rsglr/sqrt(norm(geot.ε)/ε_0)
        # 对源四面体循环
        for sid in nearCubesGeoID
            
            !(sid in geoInterval) && continue
            # 源四面体
            geos    =   geosInfolw[sid]
            # 源四面体介质对比度
            κₛ  =   tetras.κ
            # 场源距离
            Rts     =   dist(geot.center, geos.center)

            # 判断二者远近，调用不同精度的矩阵元处理函数
            if tid == sid
                # 重合
                Zts, ZtsPV  =   EFIEOnTetraPWCSepPV(geot)
                for ni in 1:3
                    n = geos.inBfsID[ni]
                    for mi in 1:3
                        # 基函数id
                        m = geot.inBfsID[mi]
                        # 写入
                        if discreteJ
                            Znear[m, n]  =   Zts[mi, ni]
                        else
                            Znear[m, n]  =   Zts[mi, ni]*κₜ
                        end
                    end
                    if discreteJ
                        Znear[n, n] += ZtsPV/(geot.ε - ε_0)
                    else
                        Znear[n, n] += ZtsPV/geot.ε
                    end
                end

            elseif Rts < Rsglrlc
                # 近场源四面体
                Zts =   EFIEOnNearTetrasPWC(geot, geos)
                # 写入数据
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geos.inBfsID[ni]
                    # 写入
                    if discreteJ
                        Znear[m, n]  =   Zts[mi, ni]
                    else
                        Znear[m, n]  =   Zts[mi, ni]*κₛ
                    end
                end
            else
                # 正常高斯求积
                # 计算四面体相关的(3*3)个矩阵元的结果
                Zts =   EFIEOnTetrasPWC(geot, geos)
                # 写入数据
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geos.inBfsID[ni]
                    # 写入
                    if discreteJ
                        Znear[m, n]  =   Zts[mi, ni]
                    else
                        Znear[m, n]  =   Zts[mi, ni]*κₛ
                    end
                end
            end # if
        end #jGeo
    end #iGeo
    return nothing
end

## TODO
"""
采用 RBF 基函数计算指定盒子内六面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在Znear稀疏矩阵中
"""
function calZnearChunkEFIEonCube!(iCube::Int, cubes, 
    geosInfo::AbstractVector{HexahedraInfo{IT, FT, CT}},
    Znear, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:LinearBasisFunction}
    
    # 本盒子信息
    cube    =   cubes[iCube]
    # 常数
    Rsglr   =   Params.Rsglr
    lockZ   =   SpinLock()

    # 盒子里的基函数区间
    cubeBFinterval  =   cube.bfInterval
    # 找出对应的三角形id
    cubeGeoID       =   cube.geoIDs

    # 找出所有邻盒子包含的基函数总数、id
    # nNearCubeBFs        =   0 
    # nearCubeBFinterval  =   Vector{UnitRange{IT}}(undef, length(cube.neighbors))
    # 邻盒子中的几何体数
    nNearCubeGeos       =   0 
    # 实际用到的邻盒子数
    nNearCubeCal    =   0
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        nNearCubeCal += 1
        # 邻盒子
        nearCube    =   cubes[jNearCube]

        # # 盒子里的基函数索引
        # nearCubeBFinterval[nNearCubeCal]   =   nearCube.bfInterval
        # # 累加基函数数量
        # nNearCubeBFs    +=  length(nearCube.bfInterval)
        # 累加计算几何体数量
        nNearCubeGeos   +=  length(nearCube.geoIDs)
        
    end
    # nearCubeBFinterval 是按所有邻盒子大小初始化的，此处要调整
    # resize!(nearCubeBFinterval, nNearCubeCal)
    # 收集所有邻盒子（不包括 id 更大的）内的基函数
    # nearCubeBFindices::Vector{IT}   =   vcat(nearCubeBFinterval...)
    # # 排序以便后面更好地使用
    # sort!(nearCubeBFindices)
    nearCubeBFindices = Znear.colIndices

    # 找出对应的三角形id
    nearCubesGeoID  =   zeros(IT, nNearCubeGeos)
    nNearCubeGeoPtr =   1
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        # 邻盒子
        nearCube    =   cubes[jNearCube]
        # 写入邻三角形
        nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs)-1)]  .=   nearCube.geoIDs
        nNearCubeGeoPtr  +=   length(nearCube.geoIDs)
    end
    # 排序并剔除冗余元素
    unique!(sort!(nearCubesGeoID))
    # fetch到本地
    geosInfolw  =   fetchDVec2LW(geosInfo, nearCubesGeoID)
    # 网格区间
    geoInterval =   getGeosInterval(geosInfo)

    # 对场盒子内四面体循环
    # @inbounds for iGeo in 1:length(nearCubesGeoID)
    @inbounds for tid in cubeGeoID
        
        # 局域的场四面体
        !(tid in geoInterval) && continue
        geot  =   geosInfolw[tid]
        #= 场四面体与源四面体在不在一个盒子？因为程序利用了目标的EFIE矩阵的对称性
        进行对称位置阻抗矩阵元的计算，要避免对同一个盒子内阻抗矩阵元的重复计算 =#
        # tins  =   nearCubesGeoID[iGeo] in cubeGeoID
        # tins  =   !isempty(searchsorted(cubeGeoID, nearCubesGeoID[iGeo]))
        # 测试四面体包含的四个测试基函数是否在当前盒子（测试盒子）的基函数（测试基函数）区间内
        msInInterval    =   [m in cubeBFinterval for m in geot.inBfsID]
        # 局部判断奇异性距离
        Rsglrlc =   Rsglr/sqrt(norm(geot.ε)/ε_0)
        # 对源四面体循环
        for sid in nearCubesGeoID
            
            !(sid in geoInterval) && continue
            # 源四面体
            geos    =   geosInfolw[sid]
            # 场源距离
            Rts     =   dist(geot.center, geos.center)
            # 源四面体包含的四个源基函数是否在所有邻盒子（源盒子）的基函数（源基函数）区间内
            nsInInterval    =   [!isempty(searchsorted(nearCubeBFindices, n)) for n in geos.inBfsID]

            # 判断二者远近，调用不同精度的矩阵元处理函数
            if tid == sid
                # 重合场源六面体
                Zts     =   EFIEOnHexaRBF(geot)
                # 写入数据

                for ni in 1:6, mi in 1:6
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geot.inBfsID[ni]
                    # 往矩阵填充结果
                    # 判断是不是在源盒子、场盒子包含的区间内
                    (msInInterval[mi] && nsInInterval[ni]) && begin
                        lock(lockZ)
                        Znear[m, n] += Zts[mi, ni]
                        unlock(lockZ)
                    end # begin
                end

            elseif Rts < Rsglrlc
                # 需要进行近奇异性处理的场源六面体
                Zts, Zst    =   EFIEOnNearHexasRBF(geot, geos)
                # 写入数据
                for ni in 1:6, mi in 1:6
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geos.inBfsID[ni]
                    # 判断是不是在源盒子、场盒子包含的区间内
                    ((msInInterval[mi] && nsInInterval[ni])) && begin
                        lock(lockZ)
                        Znear[m, n] += Zts[mi, ni]
                        unlock(lockZ)
                    end
                end
            else
                # 正常高斯求积
                # 计算六面体相关的(4*4)个矩阵元的结果
                Zts, Zst    =   EFIEOnHexasRBF(geot, geos)
                # 写入数据
                for ni in 1:6, mi in 1:6
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geos.inBfsID[ni]
                    # 判断是不是在源盒子、场盒子包含的区间内
                    ((msInInterval[mi] && nsInInterval[ni])) && begin
                        lock(lockZ)
                        Znear[m, n] += Zts[mi, ni]
                        unlock(lockZ)
                    end
                end
                
            end # if
        end #jGeo
    end #iGeo

    return nothing
end

"""
采用 PWC 基函数计算指定盒子内六面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在Znear稀疏矩阵中
"""
function calZnearChunkEFIEonCube!(iCube::Int, cubes, 
    geosInfo::AbstractVector{HexahedraInfo{IT, FT, CT}},
    Znear, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:ConstBasisFunction}
    
    # 常数
    Rsglr       =   Params.Rsglr
    # 判断体电流的离散方式
    discreteJ::Bool =   SimulationParams.discreteVar === "J"
    # 几何信息索引区间
    geoInterval =   getGeosInterval(geosInfo)

    # 本盒子信息
    cube    =   cubes[iCube]
    # 常数
    Rsglr       =   Params.Rsglr

    # 找出对应的几何体id
    cubeGeoID       =   cube.geoIDs

    # 找出所有邻盒子包含的基函数总数、id
    # nNearCubeBFs        =   0 
    # nearCubeBFinterval  =   Vector{UnitRange{IT}}(undef, length(cube.neighbors))
    # 邻盒子中的几何体数
    nNearCubeGeos       =   0 
    # 实际用到的邻盒子数
    nNearCubeCal    =   0
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        nNearCubeCal += 1
        # 邻盒子
        nearCube    =   cubes[jNearCube]
        # 累加计算几何体数量
        nNearCubeGeos   +=  length(nearCube.geoIDs)
        
    end

    # 找出对应的几何体id
    nearCubesGeoID  =   zeros(IT, nNearCubeGeos)
    nNearCubeGeoPtr =   1
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        # 邻盒子
        nearCube    =   cubes[jNearCube]
        # 写入邻几何体
        nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs)-1)]  .=   nearCube.geoIDs
        nNearCubeGeoPtr  +=   length(nearCube.geoIDs)
    end
    # 排序并剔除冗余元素
    unique!(sort!(nearCubesGeoID))
    # 将数据拉取到本地
    geosInfolw  =   fetchDVec2LW(geosInfo, nearCubesGeoID)
    
    # 对盒子内几何体循环
    @inbounds for tid in cubeGeoID
        # 局域的场几何体
        !(tid in geoInterval) && continue
        geot  =   geosInfolw[tid]
        # 场几何体介质对比度
        κₜ  =   geot.κ
        # 局部判断奇异性距离
        Rsglrlc =   Rsglr/sqrt(norm(geot.ε)/ε_0)
        # 对源几何体循环
        for sid in nearCubesGeoID
            
            !(sid in geoInterval) && continue
            # 源几何体
            geos    =   geosInfolw[sid]
            # 源几何体介质对比度
            κₛ      =   geos.κ
            # 场源距离
            Rts     =   dist(geot.center, geos.center)

            # 判断二者远近，调用不同精度的矩阵元处理函数
            if tid == sid
                # 重合
                Zts, ZtsPV    =   EFIEOnHexaPWCSepPV(geot)
                for ni in 1:3
                    # 基函数id
                    n = geos.inBfsID[ni]
                    for mi in 1:3
                        # 基函数id
                        m = geot.inBfsID[mi]
                        
                        # 写入
                        if discreteJ
                            Znear[m, n]  =   Zts[mi, ni]
                        else
                            Znear[m, n]  =   Zts[mi, ni]*κₜ
                        end
                    end
                    if discreteJ
                        Znear[n, n] += ZtsPV/(geot.ε - ε_0)
                    else
                        Znear[n, n] += ZtsPV/geot.ε
                    end
                end

            elseif Rts < Rsglrlc
                # 近场源六面体
                Zts =   EFIEOnNearHexasPWC(geot, geos)
                # 写入数据
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geos.inBfsID[ni]
                    # 写入
                    if discreteJ
                        Znear[m, n]  =   Zts[mi, ni]
                    else
                        Znear[m, n]  =   Zts[mi, ni]*κₛ
                    end
                end
            else
                # 正常高斯求积
                # 计算六面体相关的(3*3)个矩阵元的结果
                Zts =   EFIEOnHexasPWC(geot, geos)
                # 写入数据
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = geot.inBfsID[mi]
                    n = geos.inBfsID[ni]
                    # 写入
                    if discreteJ
                        Znear[m, n]  =   Zts[mi, ni]
                    else
                        Znear[m, n]  =   Zts[mi, ni]*κₛ
                    end
                end
                
            end # if
        end #jGeo
    end #iGeo

    return nothing
end


## TODO
"""
采用 PWC + PWC 基函数计算 四面体 + 六面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在Znear稀疏矩阵中
"""
function calZnearChunkEFIEonCube!(iCube::Int, cubes, 
    geos1Info::AbstractVector{VT1}, geos2Info::AbstractVector{VT2},
    Znear, ::Type{BFT}) where {VT1<:TetrahedraInfo, VT2<:HexahedraInfo, BFT<:PWC}

    # 两中网格的索引区间
    geos1Interval   =   eachindex(geos1Info)
    # geos2Interval   =   eachindex(geos2Info)
    # 离散的是否为电流
    discreteJ::Bool =   (SimulationParams.discreteVar === "J")

    # 本盒子信息
    cube    =   cubes[iCube]
    # 常数
    Rsglr       =   Params.Rsglr

    # 找出对应的三角形id
    cubeGeoID       =   cube.geoIDs

    # 找出所有邻盒子包含的基函数总数、id
    # nNearCubeBFs        =   0 
    # nearCubeBFinterval  =   Vector{UnitRange{Int}}(undef, length(cube.neighbors))
    # 邻盒子中的几何体数
    nNearCubeGeos       =   0 
    # 实际用到的邻盒子数
    nNearCubeCal    =   0
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        nNearCubeCal += 1
        # 邻盒子
        nearCube    =   cubes[jNearCube]

        # # 盒子里的基函数索引
        # nearCubeBFinterval[nNearCubeCal]   =   nearCube.bfInterval
        # # 累加基函数数量
        # nNearCubeBFs    +=  length(nearCube.bfInterval)
        # 累加计算几何体数量
        nNearCubeGeos   +=  length(nearCube.geoIDs)
        
    end
    # nearCubeBFinterval 是按所有邻盒子大小初始化的，此处要调整
    # resize!(nearCubeBFinterval, nNearCubeCal)
    # 收集所有邻盒子（不包括 id 更大的）内的基函数
    # nearCubeBFindices::Vector{Int}   =   vcat(nearCubeBFinterval...)
    # # 排序以便后面更好地使用
    # sort!(nearCubeBFindices)
    nearCubeBFindices = Znear.colIndices

    # 找出对应的三角形id
    nearCubesGeoID  =   zeros(Int, nNearCubeGeos)
    nNearCubeGeoPtr =   1
    @inbounds for j in eachindex(cube.neighbors)
        jNearCube   =   cube.neighbors[j]
        ## 分布式编程为减少数据通信，不再利用对称性填充
        # 由矩阵对称性可跳过编号较小的盒子
        # jNearCube >  iCube && continue
        # 邻盒子
        nearCube    =   cubes[jNearCube]
        # 写入邻三角形
        nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs)-1)]  .=   nearCube.geoIDs
        nNearCubeGeoPtr  +=   length(nearCube.geoIDs)
    end
    # 排序并剔除冗余元素
    unique!(sort!(nearCubesGeoID))
    # fetch 到本地
    geos1Infolw  =   fetchDVec2LW(geos1Info, filter(x -> x in geos1Interval, nearCubesGeoID))
    geos2Infolw  =   fetchDVec2LW(geos2Info, filter(x -> x in geos2Interval, nearCubesGeoID))
    # 对场盒子内三角形循环
    @inbounds for tid in cubeGeoID
                
        # 场分geo1 geo2
        tin1    =   (tid in geos1Interval)
        tin1 ? begin
            geot1    =   geos1Infolw[tid]
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(geot1.ε)/ε_0)
            # 中心
            centert =   geot1.center
        end : begin
            geot2    =   geos2Infolw[tid]
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(geot2.ε)/ε_0)
            # 中心
            centert =   geot2.center
        end
        
        # 对源体循环
        for sid in nearCubesGeoID

            # 场网格元所在的盒子为本盒子时，在 场网格元id 小于 源网格元 id时， 跳过，避免重复计算
            # tid > sid && continue
            # 源亦分geo1 geo2
            sin1    =   (sid in geos1Interval)
            sin1 ? begin
                geos1    =   geos1Infolw[sid]
                # 局部判断奇异性距离
                Rsglrlc =   Rsglr/sqrt(norm(geos1.ε)/ε_0)
                # 中心
                centers =   geos1.center
            end : begin
                geos2    =   geos2Infolw[sid]
                # 局部判断奇异性距离
                Rsglrlc =   Rsglr/sqrt(norm(geos2.ε)/ε_0)
                # 中心
                centers =   geos2.center
            end
        
            # 场源距离
            Rts     =   dist(centert, centers)
            if tid == sid
                # 重合
                Zts, ZtsPV  =   tin1 ? EFIEOnTetraPWCSepPV(geot1) : EFIEOnHexaPWCSepPV(geot2)
                tin1 ? writeZttdst!(Znear, Zts, ZtsPV, geot1, discreteJ) : writeZtt!(Znear, Zts, ZtsPV, geot2, discreteJ)

            else
                # 需要进行近奇异性处理的场源网格元
                if tin1
                    if sin1
                        Zts    =   (Rts < Rsglrlc) ? EFIEOnNearTetrasPWC(geot1, geos1) : EFIEOnTetrasPWC(geot1, geos1)
                        writeZtsdst!(Znear, Zts, geot1, geos1, discreteJ)
                    else
                        Zts    =   (Rts < Rsglrlc) ? EFIEOnNearHexaTetraPWC(geot1, geos2) : EFIEOnHexaTetraPWC(geot1, geos2)
                        writeZtsdst!(Znear, Zts, geot1, geos2, discreteJ)
                    end
                else
                    if sin1
                        Zts    =   (Rts < Rsglrlc) ? EFIEOnNearHexaTetraPWC(geot2, geos1) : EFIEOnHexaTetraPWC(geot2, geos1)
                        writeZtsdst!(Znear, Zts, geot2, geos1, discreteJ)
                    else
                        # @show 1, typeof(geot2), typeof(geos2), discreteJ
                        Zts    =   (Rts < Rsglrlc) ? EFIEOnNearHexasPWC(geot2, geos2) : EFIEOnHexasPWC(geot2, geos2)
                        writeZtsdst!(Znear, Zts, geot2, geos2, discreteJ)
                    end
                end
            end # if

        end #jGeo
    end #iGeo

    return nothing
end

"""
为适应类型变化而将写入部分单独封装
"""
function writeZttdst!(Znear, Zts, ZtsPV::T, geot::GT, discreteJ::Bool) where  {T<:Number, GT<:VolumeCellType}

    for ni in 1:3
        n = geot.inBfsID[ni]
        for mi in 1:3
            # 基函数id
            m = geot.inBfsID[mi]
            # 写入
            if discreteJ
                Znear[m, n]  =   Zts[mi, ni]
            else
                Znear[m, n]  =   Zts[mi, ni]*geot.κ
            end
        end
        if discreteJ
            Znear[n, n] += ZtsPV/(geot.ε - ε_0)
        else
            Znear[n, n] += ZtsPV/geot.ε
        end
    end

    nothing
end

function writeZtsdst!(Znear, Zts, geot::GT1, geos::GT2, discreteJ::Bool) where {GT1<:VolumeCellType, GT2<:VolumeCellType}

    # 写入数据
    for ni in 1:3, mi in 1:3
        # 基函数id
        m = geot.inBfsID[mi]
        n = geos.inBfsID[ni]
        if discreteJ
            Znear[m, n]  =   Zts[mi, ni]
        else
            Znear[m, n]  =   Zts[mi, ni]*geos.κ
        end
    end

    nothing

end


"""
采用 基函数计算指定层内 EFIE 阻抗矩阵近场元并将结果放在 ZnearChunk 中 (分布式)
"""
function calZnearChunksEFIE!(cubes, geosInfo::AbstractVector{GT},
    ZnearChunks, bfT::Type{BFT}) where {GT<:VSCellType, BFT<:BasisFunctionType}
    
    # nChunks =   length(ZnearChunks)
    # @showprogress 1 "Calculating Z Chunks" @distributed for i in 1:nChunks
    #     calZnearChunkEFIEonCube!(i, cubes, geosInfo, ZnearChunks[i], bfT)
    # end
    with_workers(calZnearChunksEFIELW!, cubes, geosInfo, ZnearChunks, bfT; withpmeter = true)

    nothing

end # function

"""
采用 基函数计算指定层内 EFIE 阻抗矩阵近场元并将结果放在 ZnearChunk 中 (分布式)
"""
function calZnearChunksEFIELW!(cubes, geosInfo::AbstractVector{GT},
    ZnearChunks, bfT::Type{BFT}) where {GT<:VSCellType, BFT<:BasisFunctionType}
    
    # 本进程索引
    idcs    =   localindices(ZnearChunks)[1]
    # 将数据 fetch 到本地
    cubeslw     =   getLocalDArgs(cubes)
    # 本进程数据
    ZnearChunkslc   =   OffsetVector(localpart(ZnearChunks), idcs[1] - 1)
    # 进度条
    cond = true#myid() == 2
    if cond
        pmeter = Progress(length(idcs); desc = "Cal Z on worker $(myid())...", dt = 1, color = :origin)
    end
    # 计算
    @threads for i in idcs
        calZnearChunkEFIEonCube!(i, cubeslw, geosInfo, ZnearChunkslc[i], bfT)
        cond && next!(pmeter)
    end

    nothing

end # function

"""
采用 基函数计算指定层内 EFIE 阻抗矩阵近场元并将结果放在 ZnearChunk 中 (分布式)
"""
function calZnearChunksEFIE!(cubes, geosInfo1::AbstractVector{T1}, geosInfo2::AbstractVector{T2},
    ZnearChunks, bfT::Type{BFT}) where {T1 <: VSCellType, T2 <: VSCellType, BFT<:BasisFunctionType}
    
    with_workers(calZnearChunksEFIELW!, cubes, geosInfo1, geosInfo2, ZnearChunks, bfT)

    nothing

end # function

"""
采用 基函数计算指定层内 EFIE 阻抗矩阵近场元并将结果放在 ZnearChunk 中 (分布式)
"""
function calZnearChunksEFIELW!(cubes, geosInfo1::AbstractVector{T1}, geosInfo2::AbstractVector{T2},
    ZnearChunks, bfT::Type{BFT}) where {T1 <: VSCellType, T2 <: VSCellType, BFT<:BasisFunctionType}
    # 本进程索引
    idcs    =   localindices(ZnearChunks)[1]
    # 将数据 fetch 到本地
    cubeslw     =   getLocalDArgs(cubes)
    # 本进程数据
    ZnearChunkslc   =   OffsetVector(localpart(ZnearChunks), idcs[1] - 1)
    # 进度条
    cond = myid() == 2
    if cond
        pmeter = Progress(length(idcs); desc = "Cal Z on worker $(myid())...", dt = 1, color = :origin)
    end
    # 计算
    @threads for i in idcs
        calZnearChunkEFIEonCube!(i, cubeslw, geosInfo1, geosInfo2, ZnearChunkslc[i], bfT)
        cond && next!(pmeter)
    end

    nothing

end # function


"""
采用 基函数计算指定层内 EFIE 阻抗矩阵近场元并将结果放在 ZnearChunk 中 (分布式)
"""
function calZnearChunksEFIET!(cubes, geosInfo::AbstractVector{GT},
    ZnearChunks, bfT::Type{BFT}) where {GT<:VSCellType, BFT<:BasisFunctionType}
    
    nChunks =   length(ZnearChunks)
    pmeter  =   Progress(nChunks, "Calculating Z Chunks")
    # 计算
    @threads for i in 1:nChunks
        calZnearChunkEFIEonCube!(i, cubes, geosInfo, ZnearChunks[i], bfT)
        next!(pmeter)
    end

    nothing

end # function

"""
采用 基函数计算指定层内 EFIE 阻抗矩阵近场元并将结果放在 ZnearChunk 中 (分布式)
"""
function calZnearChunksEFIET!(cubes, geosInfo1::AbstractVector{GT1}, geosInfo2::AbstractVector{GT2},
    ZnearChunks, bfT::Type{BFT}) where {GT1<:VSCellType, GT2<:VSCellType, BFT<:BasisFunctionType}
    
    nChunks =   length(ZnearChunks)
    pmeter  =   Progress(nChunks, "Calculating Z Chunks")
    @threads for i in 1:nChunks
        calZnearChunkEFIEonCube!(i, cubes, geosInfo1, geosInfo2, ZnearChunks[i], bfT)
        next!(pmeter)
    end

    nothing

end # function