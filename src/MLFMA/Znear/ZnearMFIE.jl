"""
采用 RWG 基函数计算 MFIE 面积分（SIE）阻抗矩阵近场元并将结果放在ZnearCSC稀疏矩阵中
"""
function calZnearCSCMFIE!(level, trianglesInfo::Vector{TriangleInfo{IT, FT}},
                        ZnearCSC::ZnearT{CT}, ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:RWG}
    
    # 本层盒子信息
    cubes   =   level.cubes
    # cubes索引
    cubesIndices    =   eachindex(cubes)
    # 常数
    Rsglr       =   Params.Rsglr
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
        @inbounds for iTri in eachindex(view(trianglesInfo, nearCubesTriID))
            # 局域的场三角形
            tid   =   nearCubesTriID[iTri]
            trit  =   trianglesInfo[tid]
            #= 场三角形与源三角形在不在一个盒子？因为程序利用了PEC目标的MFIE矩阵的对称性
            进行对称位置阻抗矩阵元的计算，要避免对同一个盒子内阻抗矩阵元的重复计算 =#
            # tins  =   nearCubesTriID[iTri] in cubeTriID
            tins  =   !isempty(searchsorted(cubeTriID, nearCubesTriID[iTri]))
            # 测试三角形包含的三个测试基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
            msInInterval    =   [!isempty(searchsorted(nearCubeBFindices, m)) for m in trit.inBfsID]

            for jTri in eachindex(view(trianglesInfo, cubeTriID))
                
                sid =   cubeTriID[jTri]
                # 源三角形
                tris    =   trianglesInfo[sid]
                # 场源距离
                Rts     =   dist(trit.center, tris.center)
                # 源三角形包含的三个源基函数是否在所有邻盒子（测试盒子）的基函数（测试基函数）区间内
                nsInInterval    =   [n in cubeBFinterval for n in tris.inBfsID]
                # 判断二者远近，调用不同精度的矩阵元处理函数
                if tid == sid
                    # 计算三角形相关的(3*3)个矩阵元的结果
                    Zts  =  MFIEOnTris(trit)
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
                    Zts, Zst    =   MFIEOnNearTris(trit, tris)
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
                            ZnearCSC[n, m] += Zst[ni, mi]
                        end
                    end
                else
                    # 正常高斯求积
                    # 计算三角形相关的(3*3)个矩阵元的结果
                    Zts, Zst    =   MFIEOnTris(trit, tris)
                    
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
                            ZnearCSC[n, m] += Zst[ni, mi]
                        end
                    end
                    
                end # if
            end #jTri
        end #iTri
    end #iCube

    return nothing
end