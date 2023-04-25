"""
采用 RWG 基函数计算 EFIE 面积分（SIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIE!(level, trianglesInfo::Vector{TriangleInfo{IT, FT}},
                        ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:LinearBasisFunction}
    
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr       =   Params.Rsglr
    ntri        =   length(trianglesInfo)
    # 叶层盒子数量
    nCubes      =   cubesIndices.stop
    # Progress Meter
    pmeter  =   Progress(nCubes; desc = "Calculating Znear (RWG)...", dt = 1)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        next!(pmeter)
        # 盒子
        cube    =   cubes[iCube]

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的三角形id
        cubeTriID       =   cube.geoIDs

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
            nNearCubeBFs    +=  length(nearCube.bfInterval)
            # 累加计算几何体数量
            nNearCubeGeos   +=  length(nearCube.geoIDs)
            
        end
        # nearCubeBFinterval 是按所有邻盒子大小初始化的，此处要调整
        resize!(nearCubeBFinterval, nNearCubeCal)
        # 收集所有邻盒子（不包括 id 更大的）内的基函数
        nearCubeBFindices::Vector{IT}   =   vcat(nearCubeBFinterval...)
        # 排序以便后面更好地使用
        sort!(nearCubeBFindices)

        # 找出对应的三角形id
        nearCubesTriID  =   zeros(IT, nNearCubeGeos)
        nNearCubeTriPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # 由矩阵对称性可跳过编号较小的盒子
            jNearCube >  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻三角形
            nearCubesTriID[nNearCubeTriPtr:(nNearCubeTriPtr+length(nearCube.geoIDs)-1)]  .=   nearCube.geoIDs
            nNearCubeTriPtr  +=   length(nearCube.geoIDs)
        end
        # 排序并剔除冗余元素
        unique!(sort!(nearCubesTriID))

        # 对盒子内三角形循环
        @inbounds for iTri in 1:length(nearCubesTriID)
            # 局域的场三角形
            tid   =   nearCubesTriID[iTri]
            # 超出边界则跳过
            tid   >   ntri && continue
            trit  =   trianglesInfo[tid]
            #= 场三角形与源三角形在不在一个盒子？因为程序利用了PEC目标的EFIE矩阵的对称性
            进行对称位置阻抗矩阵元的计算，要避免对同一个盒子内阻抗矩阵元的重复计算 =#
            # tins  =   nearCubesTriID[iTri] in cubeTriID
            # tins  =   !isempty(searchsorted(cubeTriID, nearCubesTriID[iTri]))
            # 测试三角形包含的三个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
            msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in trit.inBfsID]

            for jTri in 1:length(cubeTriID)
                
                sid =   cubeTriID[jTri]
                # 超出边界则跳过
                sid   >   ntri && continue
                # 源三角形
                tris    =   trianglesInfo[sid]
                # 场源距离
                Rts     =   dist(trit.center, tris.center)
                # 源三角形包含的三个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                nsInInterval    =   [n in cubeBFinterval for n in tris.inBfsID]
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
                            ZnearCSC[m, n] += Zts[mi, ni]
                            !(m in cubeBFinterval) && begin
                                ZnearCSC[n, m] += Zts[mi, ni]
                            end
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
                        # 避免线程锁的矩阵元循环方式下产生的条件
                        (tid > sid) && (m in cubeBFinterval) && continue

                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] & nsInInterval[ni])) && begin
                            ZnearCSC[m, n] += Zts[mi, ni]
                            ZnearCSC[n, m] += Zts[mi, ni]
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
                            ZnearCSC[m, n] += Zts[mi, ni]
                            ZnearCSC[n, m] += Zts[mi, ni]
                        end
                    end
                    
                end # if
            end #jTri
        end #iTri
    end #iCube

    return nothing
end

"""
采用 SWG 基函数计算网格元 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIE!(level, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}},
                            ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:LinearBasisFunction}
    
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr   =   Params.Rsglr
    # 是否为偏置数组（用于混合网格）
    ntetra      =   length(tetrasInfo)
    isoffset    =   isa(tetrasInfo, OffsetVector)
    geoInterval =   begin
        isoffset ? begin
            st  =   (eachindex(tetrasInfo).offset) + 1
            st:(st - 1 + ntetra)
        end : begin
            1:ntetra
        end
    end
    # 叶层盒子数量
    nCubes  =   cubesIndices.stop
    # Progress Meter
    pmeter  =   Progress(nCubes; desc = "Calculating Znear (SWG)...", dt = 1)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        next!(pmeter)
        # 盒子
        cube    =   cubes[iCube]

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的四面体id
        cubeTetraID     =   cube.geoIDs
        # cubeTetraID     =   zeros(IT, 2*length(cubeBFinterval))
        # @inbounds for (i, swg) in enumerate(view(swgsInfo, cubeBFinterval))
        #     cubeTetraID[(2i-1):(2i)]  .=  swg.inGeo
        # end
        # # 排序并剔除冗余元素
        # unique!(sort!(cubeTetraID))
        # # 防止边界元基函数的负部，编号为 0 的盒子出现在索引中
        # cubeTetraID[1] == 0 && popfirst!(cubeTetraID)
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
        
        # # 找出对应的四面体id
        # nearCubesTetraID      =   zeros(IT, 2*nNearCubeBFs)
        # @inbounds for (i, swg) in enumerate(view(swgsInfo, nearCubeBFindices))
        #     nearCubesTetraID[(2i-1):(2i)]  .=  swg.inGeo
        # end
        # 找出对应的四面体id
        nearCubesTetraID  =   zeros(IT, nNearCubeGeos)
        nNearCubeTetraPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # 由矩阵对称性可跳过编号较小的盒子
            jNearCube >  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻四面体
            nearCubesTetraID[nNearCubeTetraPtr:(nNearCubeTetraPtr+length(nearCube.geoIDs) - 1)]  .=   nearCube.geoIDs
            nNearCubeTetraPtr   +=  length(nearCube.geoIDs)
        end
        # 排序并剔除冗余元素
        unique!(sort!(nearCubesTetraID))
        # 防止边界元的编号为 0 的盒子出现在索引中
        # nearCubesTetraID[1] == 0 && popfirst!(nearCubesTetraID)
        # 对场盒子内四面体循环
        @inbounds for iTetra in 1:length(nearCubesTetraID)
            
            # 局域的场四面体
            tid     =   nearCubesTetraID[iTetra]
            !(tid in geoInterval) && continue
            tetrat  =   tetrasInfo[tid]
            #= 场四面体与源四面体在不在一个盒子？因为程序利用了目标的EFIE矩阵的对称性
            进行对称位置阻抗矩阵元的计算，要避免对同一个盒子内阻抗矩阵元的重复计算 =#
            # tins  =   nearCubesTetraID[iTetra] in cubeTetraID
            # tins  =   !isempty(searchsorted(cubeTetraID, nearCubesTetraID[iTetra]))
            # 测试四面体包含的四个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
            msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in tetrat.inBfsID]
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(tetrat.ε)/ε_0)
            # 对源四面体循环
            for jTetra in 1:length(cubeTetraID)
                
                sid =   cubeTetraID[jTetra]
                !(sid in geoInterval) && continue
                # 源四面体
                tetras    =   tetrasInfo[sid]
                # 场源距离
                Rts     =   dist(tetrat.center, tetras.center)
                # 源四面体包含的四个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                nsInInterval    =   [n in cubeBFinterval for n in tetras.inBfsID]

                # 判断二者远近，调用不同精度的矩阵元处理函数
                if tid == sid
                    # 重合场源四面体
                    Zts     =   EFIEOnTetraSWG(tetrat)
                    # 写入数据
                    for ni in 1:4, mi in 1:4
                        # 基函数id
                        m = tetrat.inBfsID[mi]
                        n = tetrat.inBfsID[ni]
                        # 往矩阵填充结果
                        # 判断是不是在源盒子、场盒子包含的区间内
                        (msInInterval[mi] && nsInInterval[ni]) && begin
                            ZnearCSC[m, n] += Zts[mi, ni]
                            !(m in cubeBFinterval) && begin
                                ZnearCSC[n, m] += Zts[mi, ni]
                            end
                        end # begin
                    end

                elseif Rts < Rsglrlc
                    # 需要进行近奇异性处理的场源四面体
                    Zts, Zst    =   EFIEOnNearTetrasSWG(tetrat, tetras)
                    # 写入数据
                    for ni in 1:4, mi in 1:4
                        # 基函数id
                        m = tetrat.inBfsID[mi]
                        n = tetras.inBfsID[ni]

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
                    Zts, Zst    =   EFIEOnTetrasSWG(tetrat, tetras)
                    # 写入数据
                    for ni in 1:4, mi in 1:4
                        # 基函数id
                        m = tetrat.inBfsID[mi]
                        n = tetras.inBfsID[ni]
                        # 避免线程锁的矩阵元循环方式下产生的条件
                        (tid > sid) && (m in cubeBFinterval) && continue

                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] && nsInInterval[ni])) && begin
                            ZnearCSC[m, n] += Zts[mi, ni]
                            ZnearCSC[n, m] += Zst[ni, mi]
                        end
                    end
                    
                end # if
            end #jTetra
        end #iTetra
    end #iCube

    return nothing
end

"""
找到 geosInfo 所在的所有 cude id
"""
function getCubeIDsWithGeos(geoInterval, cubes)
    nCubes = length(cubes)
    # 先假定全都没有
    cubeWithGeos = zeros(Bool, nCubes)

    for iCube in 1:nCubes
        # 盒子
        cube    =   cubes[iCube]
        # 找出对应的四面体id
        cubeTetraID     =   cube.geoIDs
        for jTetra in 1:length(cubeTetraID)
                
            sid =   cubeTetraID[jTetra]
            (sid in geoInterval) && begin
                cubeWithGeos[iCube] = true
                break
            end
        end
    end

    if all(cubeWithGeos)
        return collect(1:nCubes)
    else
        return (1:nCubes)[cubeWithGeos]
    end

end


"""
采用 PWC 基函数计算四面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIE!(level, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}},
                        ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:ConstBasisFunction}
    
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr   =   Params.Rsglr
    # 是否为偏置数组
    ntetra      =   length(tetrasInfo)
    isoffset    =   isa(tetrasInfo, OffsetVector)
    geoInterval =   begin
        isoffset ? begin
            st  =   (eachindex(tetrasInfo).offset) + 1
            st:(st - 1 + ntetra)
        end : begin
            1:ntetra
        end
    end

    # 判断体电流的离散方式
    discreteJ::Bool =   SimulationParams.discreteVar === "J"
    # 叶层盒子数量
    nCubes      =   cubesIndices.stop
    # Progress Meter
    pmeter  =   Progress(nCubes; desc = "Calculating Znear (PWC)...", dt = 1)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        next!(pmeter)
        # 盒子
        cube    =   cubes[iCube]
        # 找出对应的四面体id
        cubeTetraID     =   cube.geoIDs
        # 邻盒子中的几何体数
        nNearCubeGeos       =   0 
        # 实际用到的邻盒子数
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 累加计算几何体数量
            nNearCubeGeos   +=  length(nearCube.geoIDs)
        end
        
        # 找出对应的四面体id
        nearCubesTetraID  =   zeros(IT, nNearCubeGeos)
        nNearCubeTetraPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # del   由矩阵对称性可跳过编号较小的盒子
            # 消除此项，仅由四面体id筛选即可！
            # jNearCube <  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻四面体
            nearCubesTetraID[nNearCubeTetraPtr:(nNearCubeTetraPtr+length(nearCube.geoIDs) - 1)]  .=   nearCube.geoIDs
            nNearCubeTetraPtr   +=  length(nearCube.geoIDs)
        end
        # 排序方便计算
        sort!(nearCubesTetraID)
        # 防止边界元的编号为 0 的盒子出现在索引中
        # nearCubesTetraID[1] == 0 && popfirst!(nearCubesTetraID)
        # 对场盒子内四面体循环
        @inbounds for iTetra in 1:length(nearCubesTetraID)
            
            # 局域的场四面体
            tid     =   nearCubesTetraID[iTetra]
            !(tid in geoInterval) && continue
            tetrat  =   tetrasInfo[tid]
            # 场四面体介质对比度
            κₜ  =   tetrat.κ
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(tetrat.ε)/ε_0)
            # 对源四面体循环
            for jTetra in 1:length(cubeTetraID)
                
                sid =   cubeTetraID[jTetra]
                !(sid in geoInterval) && continue
                # 场四面体所在的盒子为本盒子时，在 场四面体id 小于 源四面体 id时， 跳过，避免重复计算
                tid > sid && continue
                # 源四面体
                tetras    =   tetrasInfo[sid]
                # 源四面体介质对比度
                κₛ  =   tetras.κ
                # 场源距离
                Rts     =   dist(tetrat.center, tetras.center)

                # 判断二者远近，调用不同精度的矩阵元处理函数
                if tid == sid
                    # 重合
                    Zts, ZtsPV  =   EFIEOnTetraPWCSepPV(tetrat)
                    for ni in 1:3
                        n = tetras.inBfsID[ni]
                        for mi in 1:3
                            # 基函数id
                            m = tetrat.inBfsID[mi]
                            # 写入
                            if discreteJ
                                ZnearCSC[m, n]  =   Zts[mi, ni]
                            else
                                ZnearCSC[m, n]  =   Zts[mi, ni]*κₜ
                            end
                        end
                        if discreteJ
                            ZnearCSC[n, n] += ZtsPV/(tetrat.ε - ε_0)
                        else
                            ZnearCSC[n, n] += ZtsPV/tetrat.ε
                        end
                    end

                elseif Rts < Rsglrlc
                    # 近场源四面体
                    Zts =   EFIEOnNearTetrasPWC(tetrat, tetras)
                    # 写入数据
                    for ni in 1:3, mi in 1:3
                        # 基函数id
                        m = tetrat.inBfsID[mi]
                        n = tetras.inBfsID[ni]
                        # 写入
                        if discreteJ
                            ZnearCSC[m, n]  =   Zts[mi, ni]
                            ZnearCSC[n, m]  =   Zts[mi, ni]
                        else
                            ZnearCSC[m, n]  =   Zts[mi, ni]*κₛ
                            ZnearCSC[n, m]  =   Zts[mi, ni]*κₜ
                        end
                    end
                else
                    # 正常高斯求积
                    # 计算四面体相关的(3*3)个矩阵元的结果
                    Zts =   EFIEOnTetrasPWC(tetrat, tetras)
                    # 写入数据
                    for ni in 1:3, mi in 1:3
                        # 基函数id
                        m = tetrat.inBfsID[mi]
                        n = tetras.inBfsID[ni]
                        # 写入
                        if discreteJ
                            ZnearCSC[m, n]  =   Zts[mi, ni]
                            ZnearCSC[n, m]  =   Zts[mi, ni]
                        else
                            ZnearCSC[m, n]  =   Zts[mi, ni]*κₛ
                            ZnearCSC[n, m]  =   Zts[mi, ni]*κₜ
                        end
                    end
                end # if
            end #jTetra
        end #iTetra
    end #iCube

    return nothing
end

"""
采用 RBF 基函数计算六面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIEnew!(level, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}},
                        ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:LinearBasisFunction}
    
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr   =   Params.Rsglr
    # 几何信息索引区间
    geoInterval =   getGeosInterval(hexasInfo)
    # 叶层盒子数量
    nCubes      =   cubesIndices.stop
    # Progress Meter
    pmeter  =   Progress(nCubes; desc = "Calculating Znear (RBF)...", dt = 1)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        next!(pmeter)
        # 盒子
        cube    =   cubes[iCube]

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的六面体id
        cubeHexaID     =   cube.geoIDs
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
        
        # 找出对应的六面体id
        nearCubesHexaID  =   zeros(IT, nNearCubeGeos)
        nNearCubeHexaPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # 由矩阵对称性可跳过编号较小的盒子
            jNearCube >  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻六面体
            nearCubesHexaID[nNearCubeHexaPtr:(nNearCubeHexaPtr+length(nearCube.geoIDs) - 1)]  .=   nearCube.geoIDs
            nNearCubeHexaPtr   +=  length(nearCube.geoIDs)
        end
        # 排序并剔除冗余元素
        unique!(sort!(nearCubesHexaID))
        # 防止边界元的编号为 0 的盒子出现在索引中
        # nearCubesHexaID[1] == 0 && popfirst!(nearCubesHexaID)
        # 对场盒子内六面体循环
        @inbounds for iHexa in 1:length(nearCubesHexaID)
            
            # 局域的场六面体
            tid     =   nearCubesHexaID[iHexa]
            !(tid in geoInterval) && continue
            hexat  =   hexasInfo[tid]
            # 测试六面体包含的六个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
            msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in hexat.inBfsID]
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(hexat.ε)/ε_0)
            # 对源六面体循环
            for jHexa in 1:length(cubeHexaID)
                # 源六面体 id
                sid =   cubeHexaID[jHexa]
                !(sid in geoInterval) && continue
                # 源六面体
                hexas    =   hexasInfo[sid]
                # 场源距离
                Rts     =   dist(hexat.center, hexas.center)
                # 源六面体包含的六个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                nsInInterval    =   [n in cubeBFinterval for n in hexas.inBfsID]
                # 判断二者远近，调用不同精度的矩阵元处理函数
                if tid == sid
                    # 重合场源六面体
                    Zts     =   EFIEOnHexaRBF(hexat)
                    # 写入数据

                    for ni in 1:6, mi in 1:6
                        # 基函数id
                        m = hexat.inBfsID[mi]
                        n = hexat.inBfsID[ni]
                        # 往矩阵填充结果
                        # 判断是不是在源盒子、场盒子包含的区间内
                        (msInInterval[mi] && nsInInterval[ni]) && begin
                            ZnearCSC[m, n] += Zts[mi, ni]
                            !(m in cubeBFinterval) && begin
                                ZnearCSC[n, m] += Zts[mi, ni]
                            end
                        end # begin
                    end

                elseif Rts < Rsglrlc
                    # 需要进行近奇异性处理的场源六面体
                    Zts, Zst    =   EFIEOnNearHexasRBF(hexat, hexas)
                    # 写入数据
                    for ni in 1:6, mi in 1:6
                        # 基函数id
                        m = hexat.inBfsID[mi]
                        n = hexas.inBfsID[ni]
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
                    # 计算六面体相关的(4*4)个矩阵元的结果
                    Zts, Zst    =   EFIEOnHexasRBF(hexat, hexas)
                    # 写入数据
                    for ni in 1:6, mi in 1:6
                        # 基函数id
                        m = hexat.inBfsID[mi]
                        n = hexas.inBfsID[ni]
                        # 避免线程锁的矩阵元循环方式下产生的条件
                        (tid > sid) && (m in cubeBFinterval) && continue

                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] && nsInInterval[ni])) && begin
                            ZnearCSC[m, n] += Zts[mi, ni]
                            ZnearCSC[n, m] += Zst[ni, mi]
                        end
                    end
                    
                end # if
            end #jHexa
        end #iHexa
    end #iCube
    return nothing
end

"""
采用 RBF 基函数计算六面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIE!(level, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}},
                        ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:LinearBasisFunction}
    
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr   =   Params.Rsglr
    # 几何信息索引区间
    isoffset    =   isa(hexasInfo, OffsetVector)
    geoInterval =   getGeosInterval(hexasInfo)
    # 叶层盒子数量
    nCubes      =   cubesIndices.stop
    # 线程锁防止对同一数据写入出错
    lockZ   =   SpinLock()
    # Progress Meter
    pmeter  =   Progress(nCubes; desc = "Calculating Znear (RBF)...", dt = 1)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        next!(pmeter)
        # 盒子
        cube    =   cubes[iCube]

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的六面体id
        cubeHexaID     =   cube.geoIDs
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
        nearCubesHexaID  =   zeros(IT, nNearCubeGeos)
        nNearCubeHexaPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # del   由矩阵对称性可跳过编号较小的盒子
            # 消除此项，仅由六面体id筛选即可！
            # jNearCube <  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻六面体
            nearCubesHexaID[nNearCubeHexaPtr:(nNearCubeHexaPtr+length(nearCube.geoIDs) - 1)]  .=   nearCube.geoIDs
            nNearCubeHexaPtr   +=  length(nearCube.geoIDs)
        end
        # 排序并剔除冗余元素
        unique!(sort!(nearCubesHexaID))
        # 防止边界元的编号为 0 的盒子出现在索引中
        # nearCubesHexaID[1] == 0 && popfirst!(nearCubesHexaID)
        # 对场盒子内六面体循环
        @inbounds for iHexa in 1:length(nearCubesHexaID)
            
            # 局域的场六面体
            tid     =   nearCubesHexaID[iHexa]
            isoffset && begin
                !(tid in geoInterval) && continue
            end
            hexat  =   hexasInfo[tid]
            # 测试六面体包含的六个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
            msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in hexat.inBfsID]
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(hexat.ε)/ε_0)
            # 对源六面体循环
            for jHexa in 1:length(cubeHexaID)
                # 源六面体 id
                sid =   cubeHexaID[jHexa]
                isoffset && begin
                    !(sid in geoInterval) && continue
                end
                # 场六面体所在的盒子为本盒子时，在 场六面体id 小于 源六面体 id时， 跳过，避免重复计算
                tid > sid && continue
                # 源六面体
                hexas    =   hexasInfo[sid]
                # 场源距离
                Rts     =   dist(hexat.center, hexas.center)
                # 源六面体包含的六个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                nsInInterval    =   [n in cubeBFinterval for n in hexas.inBfsID]
                # 判断二者远近，调用不同精度的矩阵元处理函数
                if tid == sid
                    # 重合场源六面体
                    Zts     =   EFIEOnHexaRBF(hexat)
                    # 写入数据

                    for ni in 1:6, mi in 1:6
                        # 基函数id
                        m = hexat.inBfsID[mi]
                        n = hexat.inBfsID[ni]
                        # 往矩阵填充结果
                        # 判断是不是在源盒子、场盒子包含的区间内
                        (msInInterval[mi] && nsInInterval[ni]) && begin
                            lock(lockZ)
                            ZnearCSC[m, n] += Zts[mi, ni]
                            unlock(lockZ)
                        end # begin
                    end

                elseif Rts < Rsglrlc
                    # 需要进行近奇异性处理的场源六面体
                    Zts, Zst    =   EFIEOnNearHexasRBF(hexat, hexas)
                    # 写入数据
                    for ni in 1:6, mi in 1:6
                        # 基函数id
                        m = hexat.inBfsID[mi]
                        n = hexas.inBfsID[ni]
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
                    Zts, Zst    =   EFIEOnHexasRBF(hexat, hexas)
                    # 写入数据
                    for ni in 1:6, mi in 1:6
                        # 基函数id
                        m = hexat.inBfsID[mi]
                        n = hexas.inBfsID[ni]
                        # 判断是不是在源盒子、场盒子包含的区间内
                        ((msInInterval[mi] && nsInInterval[ni])) && begin
                            lock(lockZ)
                            ZnearCSC[m, n] += Zts[mi, ni]
                            ZnearCSC[n, m] += Zst[ni, mi]
                            unlock(lockZ)
                        end
                    end
                    
                end # if
            end #jHexa
        end #iHexa
    end #iCube

    return nothing
end

"""
采用 PWC 基函数计算六面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIE!(level, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}},
                        ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:ConstBasisFunction}
    
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr       =   Params.Rsglr
    # 判断体电流的离散方式
    discreteJ::Bool = SimulationParams.discreteVar == "J"
    # 几何信息索引区间
    geoInterval =   getGeosInterval(hexasInfo)
    # 叶层盒子数量
    nCubes      =   cubesIndices.stop
    # Progress Meter
    pmeter  =   Progress(nCubes; desc = "Calculating Znear (PWC)...", dt = 1)
    # 对盒子循环计算
    @threads for iCube in 1:nCubes
        next!(pmeter)
        # 盒子
        cube    =   cubes[iCube]
        # 找出对应的六面体id
        cubeHexaID     =   cube.geoIDs
        # 邻盒子中的几何体数
        nNearCubeGeos       =   0 
        # 实际用到的邻盒子数
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 累加计算几何体数量
            nNearCubeGeos   +=  length(nearCube.geoIDs)
        end
        
        # 找出对应的六面体id
        nearCubesHexaID  =   zeros(IT, nNearCubeGeos)
        nNearCubeHexaPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # del   由矩阵对称性可跳过编号较小的盒子
            # 消除此项，仅由六面体id筛选即可！
            # jNearCube <  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻六面体
            nearCubesHexaID[nNearCubeHexaPtr:(nNearCubeHexaPtr+length(nearCube.geoIDs) - 1)]  .=   nearCube.geoIDs
            nNearCubeHexaPtr   +=  length(nearCube.geoIDs)
        end
        # 排序方便计算
        sort!(nearCubesHexaID)
        # 防止边界元的编号为 0 的盒子出现在索引中
        # nearCubesHexaID[1] == 0 && popfirst!(nearCubesHexaID)
        # 对场盒子内六面体循环
        @inbounds for iHexa in 1:length(nearCubesHexaID)
            
            # 局域的场六面体
            tid     =   nearCubesHexaID[iHexa]
            !(tid in geoInterval) && continue
            hexat  =   hexasInfo[tid]
            # 场六面体介质对比度
            κₜ  =   hexat.κ
            # 局部判断奇异性距离
            Rsglrlc =   Rsglr/sqrt(norm(hexat.ε)/ε_0)
            # 对源六面体循环
            for jHexa in 1:length(cubeHexaID)
                # 源六面体 id
                sid =   cubeHexaID[jHexa]
                !(sid in geoInterval) && continue
                # 场六面体所在的盒子为本盒子时，在 场六面体id 小于 源六面体 id时， 跳过，避免重复计算
                tid < sid && continue
                # 源六面体
                hexas    =   hexasInfo[sid]
                # 源六面体介质对比度
                κₛ  =   hexas.κ
                # 场源距离
                Rts     =   dist(hexat.center, hexas.center)

                # 判断二者远近，调用不同精度的矩阵元处理函数
                if tid == sid
                    # 重合
                    Zts, ZtsPV    =   EFIEOnHexaPWCSepPV(hexat)
                    for ni in 1:3
                        # 基函数id
                        n = hexas.inBfsID[ni]
                        for mi in 1:3
                            # 基函数id
                            m = hexat.inBfsID[mi]
                            
                            # 写入
                            if discreteJ
                                ZnearCSC[m, n]  =   Zts[mi, ni]
                            else
                                ZnearCSC[m, n]  =   Zts[mi, ni]*κₜ
                            end
                        end
                        if discreteJ
                            ZnearCSC[n, n] += ZtsPV/(hexat.ε - ε_0)
                        else
                            ZnearCSC[n, n] += ZtsPV/hexat.ε
                        end
                    end

                elseif Rts < Rsglrlc
                    # 近场源六面体
                    Zts =   EFIEOnNearHexasPWC(hexat, hexas)
                    # 写入数据
                    for ni in 1:3, mi in 1:3
                        # 基函数id
                        m = hexat.inBfsID[mi]
                        n = hexas.inBfsID[ni]
                        # 写入
                        if discreteJ
                            ZnearCSC[m, n]  =   Zts[mi, ni]
                            ZnearCSC[n, m]  =   Zts[mi, ni]
                        else
                            ZnearCSC[m, n]  =   Zts[mi, ni]*κₛ
                            ZnearCSC[n, m]  =   Zts[mi, ni]*κₜ
                        end
                    end
                else
                    # 正常高斯求积
                    # 计算六面体相关的(3*3)个矩阵元的结果
                    Zts =   EFIEOnHexasPWC(hexat, hexas)
                    # 写入数据
                    for ni in 1:3, mi in 1:3
                        # 基函数id
                        m = hexat.inBfsID[mi]
                        n = hexas.inBfsID[ni]
                        # 写入
                        if discreteJ
                            ZnearCSC[m, n]  =   Zts[mi, ni]
                            ZnearCSC[n, m]  =   Zts[mi, ni]
                        else
                            ZnearCSC[m, n]  =   Zts[mi, ni]*κₛ
                            ZnearCSC[n, m]  =   Zts[mi, ni]*κₜ
                        end
                    end
                    
                end # if
            end #jHexa
        end #iHexa
    end #iCube

    return nothing
end

"""
采用 PWC + PWC 基函数计算 四面体 + 六面体 EFIE 的体积分（VIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCEFIE!(  level, geos1Info::AbstractVector{VT1}, geos2Info::AbstractVector{VT2},
                            ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {FT<:Real, CT<:Complex{FT}, 
                            VT1<:TetrahedraInfo, VT2<:HexahedraInfo, BFT<:PWC}

    # 两中网格的索引区间
    geos1Interval   =   eachindex(geos1Info)
    geos2Interval   =   eachindex(geos2Info)
    # 离散的是否为电流
    discreteJ::Bool =   (SimulationParams.discreteVar === "J")

    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr   =   Params.Rsglr
    # 叶层盒子数量
    nCubes  =   cubesIndices.stop
    # Progress Meter
    pmeter  =   Progress(nCubes; desc = "Calculating Znear (PWC + PWC)...", dt = 1)
    # 对盒子循环计算@threads
    @floop WorkStealingEx(; basesize = 1 ) for iCube in 1:nCubes
        next!(pmeter)
        # 盒子
        cube    =   cubes[iCube]

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的四面体id
        cubeGeoID       =   cube.geoIDs
        # 找出所有邻盒子包含的基函数总数、id
        nNearCubeBFs        =   0 
        nearCubeBFinterval  =   Vector{UnitRange{Int}}(undef, length(cube.neighbors))
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
        nearCubeBFindices::Vector{Int}   =   vcat(nearCubeBFinterval...)
        # 排序以便后面更好地使用
        sort!(nearCubeBFindices)
        
        # 找出对应的网格元id
        nearCubesGeoID  =   zeros(Int, nNearCubeGeos)
        nNearCubeGeoPtr =   1
        @inbounds for j in eachindex(cube.neighbors)
            jNearCube   =   cube.neighbors[j]
            # del   由矩阵对称性可跳过编号较小的盒子
            # 消除此项，仅由网格元id筛选即可！
            # jNearCube <  iCube && continue
            # 邻盒子
            nearCube    =   cubes[jNearCube]
            # 写入邻网格元
            nearCubesGeoID[nNearCubeGeoPtr:(nNearCubeGeoPtr+length(nearCube.geoIDs) - 1)]  .=   nearCube.geoIDs
            nNearCubeGeoPtr   +=  length(nearCube.geoIDs)
        end
        # 排序并剔除冗余元素
        unique!(sort!(nearCubesGeoID))
        # 对场盒子内网格元循环
        @inbounds for tid in nearCubesGeoID
            
            # 局域的场网格元
            # tid     =   nearCubesGeoID[iGeo]

            # 场分geo1 geo2
            tin1    =   (tid in geos1Interval)
            tin1 ? begin
                geot1    =   geos1Info[tid]
                # 测试网格元包含的四个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in geot1.inBfsID]
                # 局部判断奇异性距离
                Rsglrlc =   Rsglr/sqrt(norm(geot1.ε)/ε_0)
                # 中心
                centert =   geot1.center
            end : begin
                geot2    =   geos2Info[tid]
                # 测试网格元包含的四个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in geot2.inBfsID]
                # 局部判断奇异性距离
                Rsglrlc =   Rsglr/sqrt(norm(geot2.ε)/ε_0)
                # 中心
                centert =   geot2.center
            end
            
            # 对源体循环
            for jGeo in 1:length(cubeGeoID)
                
                sid =   cubeGeoID[jGeo]
                # 场网格元所在的盒子为本盒子时，在 场网格元id 小于 源网格元 id时， 跳过，避免重复计算
                tid > sid && continue
                # 源亦分geo1 geo2
                sin1    =   (sid in geos1Interval)
                sin1 ? begin
                    geos1    =   geos1Info[sid]
                    # 源网格元包含的四个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                    nsInInterval    =   [n in cubeBFinterval for n in geos1.inBfsID]
                    # 局部判断奇异性距离
                    Rsglrlc =   Rsglr/sqrt(norm(geos1.ε)/ε_0)
                    # 中心
                    centers =   geos1.center
                end : begin
                    geos2    =   geos2Info[sid]
                    # 源四面体包含的四个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                    nsInInterval    =   [n in cubeBFinterval for n in geos2.inBfsID]
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
                    tin1 ? writeZtt!(ZnearCSC, Zts, ZtsPV, geot1, discreteJ) : writeZtt!(ZnearCSC, Zts, ZtsPV, geot2, discreteJ)

                else
                    # 需要进行近奇异性处理的场源网格元
                    if tin1
                        if sin1
                            Zts    =   (Rts < Rsglrlc) ? EFIEOnNearTetrasPWC(geot1, geos1) : EFIEOnTetrasPWC(geot1, geos1)
                            writeZts!(ZnearCSC, Zts, msInInterval, nsInInterval, geot1, geos1, discreteJ)
                        else
                            Zts    =   (Rts < Rsglrlc) ? EFIEOnNearHexaTetraPWC(geot1, geos2) : EFIEOnHexaTetraPWC(geot1, geos2)
                            writeZts!(ZnearCSC, Zts, msInInterval, nsInInterval, geot1, geos2, discreteJ)
                        end
                    else
                        if sin1
                            Zts    =   (Rts < Rsglrlc) ? EFIEOnNearHexaTetraPWC(geot2, geos1) : EFIEOnHexaTetraPWC(geot2, geos1)
                            writeZts!(ZnearCSC, Zts, msInInterval, nsInInterval, geot2, geos1, discreteJ)
                        else
                            # @show 1, msInInterval, nsInInterval, typeof(geot2), typeof(geos2), discreteJ
                            Zts    =   (Rts < Rsglrlc) ? EFIEOnNearHexasPWC(geot2, geos2) : EFIEOnHexasPWC(geot2, geos2)
                            writeZts!(ZnearCSC, Zts, msInInterval, nsInInterval, geot2, geos2, discreteJ)
                        end
                    end
                end # if

            end #jGeo
        end #iGeo
    end #iCube

    return nothing
end

"""
为适应类型变化而将写入部分单独封装
"""
function writeZtt!(ZnearCSC, Zts, ZtsPV::T, geot::GT, discreteJ::Bool) where  {T<:Number, GT<:VolumeCellType}

    for ni in 1:3
        n = geot.inBfsID[ni]
        for mi in 1:3
            # 基函数id
            m = geot.inBfsID[mi]
            # 写入
            if discreteJ
                ZnearCSC[m, n]  =   Zts[mi, ni]
            else
                ZnearCSC[m, n]  =   Zts[mi, ni]*geot.κ
            end
        end
        if discreteJ
            ZnearCSC[n, n] += ZtsPV/(geot.ε - ε_0)
        else
            ZnearCSC[n, n] += ZtsPV/geot.ε
        end
    end

    nothing
end

function writeZts!(ZnearCSC, Zts, msInInterval, nsInInterval, geot::GT1, geos::GT2, discreteJ::Bool) where {GT1<:VolumeCellType, GT2<:VolumeCellType}

    # 写入数据
    for ni in 1:3, mi in 1:3
        # 基函数id
        m = geot.inBfsID[mi]
        n = geos.inBfsID[ni]
        # 判断是不是在源盒子、场盒子包含的区间内
        ((msInInterval[mi] && nsInInterval[ni])) && begin
            if discreteJ
                ZnearCSC[m, n]  =   Zts[mi, ni]
                ZnearCSC[n, m]  =   Zts[mi, ni]
            else
                ZnearCSC[m, n]  =   Zts[mi, ni]*geos.κ
                ZnearCSC[n, m]  =   Zts[mi, ni]*geot.κ
            end
        end
    end

    nothing

end