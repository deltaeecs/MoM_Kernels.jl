
"""
采用 RWG + SWG 基函数计算 三角形 + 四面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIE!(level, tris::Vector{TriangleInfo{IT, FT}}, tetras::AbstractVector{TetrahedraInfo{IT, FT, CT}},
                            ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:SWG}
    
    # 三角形、四面体数
    ntri    =   length(tris)
    ntetra  =   length(tetras)

    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr   =   Params.Rsglr
    # 是否为偏置数组（用于混合网格）
    ntetra      =   length(tetras)
    isoffset    =   isa(tetras, OffsetVector)
    geoInterval =   begin
        isoffset ? begin
            st  =   (eachindex(tetras).offset) + 1
            st:(st - 1 + ntetra)
        end : begin
            1:ntetra
        end
    end
    # 叶层盒子数量
    nCubes  =   cubesIndices.stop
    # Progress Meter
    pmeter  =   Progress(nCubes; desc = "Calculating Znear (RWG + SWG)...", dt = 1)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        next!(pmeter)
        # 盒子
        cube    =   cubes[iCube]

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的四面体id
        cubeGeoID       =   cube.geoIDs
        # 找出所有邻盒子包含的基函数总数、id
        nNearCubeBFs        =   0 
        nearCubeBFinterval  =   Vector{UnitRange{IT}}(undef, length(cube.neighbors))
        # 邻盒子中的几何体数
        nNearCubeGeos       =   0 
        # 实际用到的邻盒子数
        nNearCubeCal    =   0
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # 由矩阵对称性可跳过编号较小的盒子
            jNearCube >  iCube && continue
            nNearCubeCal += 1
            # 邻盒子
            nearCube    =   cubes[jNearCube]

            # 盒子里的基函数索引
            nearCubeBFinterval[nNearCubeCal]   =   nearCube.bfInterval
            # 累加基函数数量
            nNearCubeBFs    +=    length(nearCube.bfInterval)
            # 累加计算几何体数量
            nNearCubeGeos   +=  length(nearCube.geoIDs)
            
        end
        # nearCubeBFinterval 是按所有邻盒子大小初始化的，此处要调整
        resize!(nearCubeBFinterval, nNearCubeCal)
        # 收集所有邻盒子（不包括 id 更大的）内的基函数
        nearCubeBFindices::Vector{IT}   =   vcat(nearCubeBFinterval...)
        # 排序以便后面更好地使用
        sort!(nearCubeBFindices)
        
        # 找出对应的四面体id
        nearCubesGeoID  =   zeros(IT, nNearCubeGeos)
        nNearCubeGeoPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # 由矩阵对称性可跳过编号较小的盒子
            jNearCube >  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻四面体
            nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs) - 1)]  .=   nearCube.geoIDs
            nNearCubeGeoPtr   +=  length(nearCube.geoIDs)
        end
        # 排序并剔除冗余元素
        unique!(sort!(nearCubesGeoID))
        # 对场盒子内三角形循环
        @inbounds for iGeo in 1:length(nearCubesGeoID)
            
            # 局域的场四面体
            tid     =   nearCubesGeoID[iGeo]
            # 场只找三角形
            tid     >   ntri && continue
            geot    =   tris[tid]
            # 测试四面体包含的四个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
            msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in geot.inBfsID]
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(geot.ε)/ε_0)
            # 对源四面体循环
            for jGeo in 1:length(cubeGeoID)
                
                sid =   cubeGeoID[jGeo]
                # 源只找四面体
                !(sid in geoInterval) && continue
                # 源四面体
                geos    =   tetras[sid]
                # 场源距离
                Rts     =   dist(geot.center, geos.center)
                # 源四面体包含的四个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                nsInInterval    =   [n in cubeBFinterval for n in geos.inBfsID]

                if Rts < Rsglrlc
                    # 需要进行近奇异性处理的场源四面体
                    Zts, Zst    =   EFIEOnNearRWGSWG(geot, geos)
                    # 写入数据
                    for ni in 1:4, mi in 1:3
                        # 基函数id
                        m = geot.inBfsID[mi]
                        n = geos.inBfsID[ni]
                        # RWG基函数不设半基函数因此跳过
                        (m == 0) && continue
                        # 避免线程锁的矩阵元循环方式下产生的条件
                        (tid > sid) && (m in cubeBFinterval) && continue

                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] && nsInInterval[ni])) && begin
                            ZnearCSC[m, n] += Zts[mi, ni]
                            ZnearCSC[n, m] += Zst[ni, mi]
                        end
                    end
                else
                    # 正常高斯求积
                    # 计算四面体相关的(4*4)个矩阵元的结果
                    Zts, Zst    =   EFIEOnRWGSWG(geot, geos)
                    # 写入数据
                    for ni in 1:4, mi in 1:3
                        # 基函数id
                        m = geot.inBfsID[mi]
                        n = geos.inBfsID[ni]
                        # RWG基函数不设半基函数因此跳过
                        (m == 0) && continue
                        # 避免线程锁的矩阵元循环方式下产生的条件
                        (tid > sid) && (m in cubeBFinterval) && continue

                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] && nsInInterval[ni])) && begin
                            ZnearCSC[m, n] += Zts[mi, ni]
                            ZnearCSC[n, m] += Zst[ni, mi]
                        end
                    end
                    
                end # if

            end #jGeo
        end #iGeo
    end #iCube

    return nothing
end

"""
采用 RWG + RBF 基函数计算六面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIE!(level, tris::Vector{TriangleInfo{IT, FT}}, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}},
    ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:RBF}
    
    # 三角形、六面体数
    ntri    =   length(tris)
    nhexa   =   length(hexasInfo)

    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr   =   Params.Rsglr
    # 几何信息索引区间
    geoInterval =   getGeosInterval(hexasInfo)
    # 叶层盒子数量
    nCubes  =   cubesIndices.stop
    # 线程锁防止对同一数据写入出错
    lockZ   =   SpinLock()
    # Progress Meter
    Progress(nCubes; desc = "Calculating Znear (RWG + RBF)...", dt = 1)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        # 盒子
        cube    =   cubes[iCube]
        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的六面体id
        cubeGeoID       =   cube.geoIDs
        # 找出所有邻盒子包含的基函数总数、id
        nNearCubeBFs        =   0 
        nearCubeBFinterval  =   Vector{UnitRange{IT}}(undef, length(cube.neighbors))
        # 邻盒子中的几何体数
        nNearCubeGeos       =   0 
        # 实际用到的邻盒子数
        nNearCubeCal    =   0
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # del   由矩阵对称性可跳过编号较小的盒子
            # 消除此项，仅由六面体id筛选即可！
            # jNearCube <  iCube && continue
            nNearCubeCal += 1
            # 邻盒子
            nearCube    =   cubes[jNearCube]

            # 盒子里的基函数索引
            nearCubeBFinterval[nNearCubeCal]   =   nearCube.bfInterval
            # 累加基函数数量
            nNearCubeBFs    +=    length(nearCube.bfInterval)
            # 累加计算几何体数量
            nNearCubeGeos   +=  length(nearCube.geoIDs)
            
        end
        # 所有的邻盒子的基函数
        nearCubeBFindices::Vector{IT}   =   vcat(nearCubeBFinterval...)
        # 排序以便后面更好地使用
        sort!(nearCubeBFindices)
        
        # 找出对应的六面体id
        nearCubesGeoID  =   zeros(IT, nNearCubeGeos)
        nNearCubeGeoPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # del   由矩阵对称性可跳过编号较小的盒子
            # 消除此项，仅由六面体id筛选即可！
            # jNearCube <  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻六面体
            nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs) - 1)]  .=   nearCube.geoIDs
            nNearCubeGeoPtr   +=  length(nearCube.geoIDs)
        end
        # 排序并剔除冗余元素
        unique!(sort!(nearCubesGeoID))
        # 对场盒子内三角形循环
        @inbounds for iGeo in 1:length(nearCubesGeoID)
            
            # 局域的场六面体
            tid     =   nearCubesGeoID[iGeo]
            # 场只找三角形
            tid     >   ntri && continue
            geot    =   tris[tid]
            # 测试六面体包含的六个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
            msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in geot.inBfsID]
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(geot.ε)/ε_0)
            # 对源六面体循环
            for jGeo in 1:length(cubeGeoID)
                
                sid =   cubeGeoID[jGeo]
                # 源只找六面体
                !(sid in geoInterval) && continue
                # 源六面体
                geos    =   hexasInfo[sid]
                # 场源距离
                Rts     =   dist(geot.center, geos.center)
                # 源六面体包含的六个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                nsInInterval    =   [n in cubeBFinterval for n in geos.inBfsID]

                if Rts < Rsglrlc
                    # 需要进行近奇异性处理的场源六面体
                    Zts, Zst    =   EFIEOnNearRWGRBF(geot, geos)
                    # 写入数据
                    for ni in 1:6, mi in 1:3
                        # 基函数id
                        m = geot.inBfsID[mi]
                        n = geos.inBfsID[ni]
                        # RWG基函数不设半基函数因此跳过
                        (m == 0) && continue
                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] && nsInInterval[ni])) && begin
                            lock(lockZ)
                            ZnearCSC[m, n] += Zts[mi, ni]
                            ZnearCSC[n, m] += Zst[ni, mi]
                            unlock(lockZ)
                        end
                    end
                else
                    # 正常高斯求积
                    # 计算六面体相关的(4*4)个矩阵元的结果
                    Zts, Zst    =   EFIEOnRWGRBF(geot, geos)
                    # 写入数据
                    for ni in 1:6, mi in 1:3
                        # 基函数id
                        m = geot.inBfsID[mi]
                        n = geos.inBfsID[ni]
                        # RWG基函数不设半基函数因此跳过
                        (m == 0) && continue
                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] && nsInInterval[ni])) && begin
                            lock(lockZ)
                            ZnearCSC[m, n] += Zts[mi, ni]
                            ZnearCSC[n, m] += Zst[ni, mi]
                            unlock(lockZ)
                        end
                    end
                    
                end # if

            end #jGeo
        end #iGeo
        next!(pmeter)
    end #iCube

    return nothing
end

"""
采用 RWG + PWC 基函数计算 三角形 + 四面体/六面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIE!(level, tris::Vector{TriangleInfo{IT, FT}}, geosInfo::AbstractVector{VT},
                            ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, VT<:VolumeCellType, BFT<:PWC}
    
    # 三角形、四面体数
    ntri    =   length(tris)
    ngoeV  =   length(geosInfo)

    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr   =   Params.Rsglr
    # 是否为偏置数组（用于混合网格）
    isoffset    =   isa(geosInfo, OffsetVector)
    geoInterval =   begin
        isoffset ? begin
            st  =   (eachindex(geosInfo).offset) + 1
            st:(st - 1 + ngoeV)
        end : begin
            1:ngoeV
        end
    end
    # 叶层盒子数量
    nCubes  =   cubesIndices.stop
    # 线程锁防止对同一数据写入出错
    lockZ   =   SpinLock()
    # Progress Meter
    pmeter  =   Progress(nCubes; desc = "Calculating Znear (RWG + PWC)...", dt = 1)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        next!(pmeter)
        # 盒子
        cube    =   cubes[iCube]

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的四面体id
        cubeGeoID       =   cube.geoIDs
        # 找出所有邻盒子包含的基函数总数、id
        nNearCubeBFs        =   0 
        nearCubeBFinterval  =   Vector{UnitRange{IT}}(undef, length(cube.neighbors))
        # 邻盒子中的几何体数
        nNearCubeGeos       =   0 
        # 实际用到的邻盒子数
        nNearCubeCal    =   0
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # del   由矩阵对称性可跳过编号较小的盒子
            # 消除此项，仅由四面体id筛选即可！
            # jNearCube <  iCube && continue
            nNearCubeCal += 1
            # 邻盒子
            nearCube    =   cubes[jNearCube]

            # 盒子里的基函数索引
            nearCubeBFinterval[nNearCubeCal]   =   nearCube.bfInterval
            # 累加基函数数量
            nNearCubeBFs    +=    length(nearCube.bfInterval)
            # 累加计算几何体数量
            nNearCubeGeos   +=  length(nearCube.geoIDs)
            
        end
        # 所有的邻盒子的基函数
        nearCubeBFindices::Vector{IT}   =   vcat(nearCubeBFinterval...)
        # 排序以便后面更好地使用
        sort!(nearCubeBFindices)
        
        # 找出对应的四面体id
        nearCubesGeoID  =   zeros(IT, nNearCubeGeos)
        nNearCubeGeoPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # del   由矩阵对称性可跳过编号较小的盒子
            # 消除此项，仅由四面体id筛选即可！
            # jNearCube <  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻四面体
            nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs) - 1)]  .=   nearCube.geoIDs
            nNearCubeGeoPtr   +=  length(nearCube.geoIDs)
        end
        # 排序并剔除冗余元素
        unique!(sort!(nearCubesGeoID))
        # 对场盒子内三角形循环
        @inbounds for iGeo in 1:length(nearCubesGeoID)
            
            # 局域的场四面体
            tid     =   nearCubesGeoID[iGeo]
            # 场只找三角形
            tid     >   ntri && continue
            geot    =   tris[tid]
            # 测试四面体包含的四个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
            msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in geot.inBfsID]
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(geot.ε)/ε_0)
            # 对源四面体循环
            for jGeo in 1:length(cubeGeoID)
                
                sid =   cubeGeoID[jGeo]
                # 源只找四面体
                !(sid in geoInterval) && continue
                # 源四面体
                geos    =   geosInfo[sid]
                κs      =   geos.κ
                # 场源距离
                Rts     =   dist(geot.center, geos.center)
                # 源四面体包含的四个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                nsInInterval    =   [n in cubeBFinterval for n in geos.inBfsID]

                if Rts < Rsglrlc
                    # 需要进行近奇异性处理的场源四面体
                    Zts, Zst    =   EFIEOnNearRWGPWC(geot, geos)
                    # 写入数据
                    for ni in 1:3, mi in 1:3
                        # 基函数id
                        m = geot.inBfsID[mi]
                        n = geos.inBfsID[ni]
                        # RWG基函数不设半基函数因此跳过
                        (m == 0) && continue
                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] && nsInInterval[ni])) && begin
                            lock(lockZ)
                            if discreteJ
                                ZnearCSC[m, n]  +=  Zts[mi, ni]
                            else
                                ZnearCSC[m, n]  +=  Zts[mi, ni] * κs
                            end
                            ZnearCSC[n, m] += Zst[ni, mi]
                            unlock(lockZ)
                        end
                    end
                else
                    # 正常高斯求积
                    # 计算四面体相关的(4*4)个矩阵元的结果
                    Zts, Zst    =   EFIEOnRWGPWC(geot, geos)
                    # 写入数据
                    for ni in 1:3, mi in 1:3
                        # 基函数id
                        m = geot.inBfsID[mi]
                        n = geos.inBfsID[ni]
                        # RWG基函数不设半基函数因此跳过
                        (m == 0) && continue
                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] && nsInInterval[ni])) && begin
                            lock(lockZ)
                            if discreteJ
                                ZnearCSC[m, n]  +=  Zts[mi, ni]
                            else
                                ZnearCSC[m, n]  +=  Zts[mi, ni] * κs
                            end
                            ZnearCSC[n, m] += Zst[ni, mi]
                            unlock(lockZ)
                        end
                    end
                    
                end # if

            end #jGeo
        end #iGeo
    end #iCube

    return nothing
end