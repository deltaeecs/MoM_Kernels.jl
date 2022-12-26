
"""
在叶层从基函数向盒子聚合
"""
function aggOnBF!(level, aggSBF, IVec::AbstractArray{T}) where T<:Number
    # 叶层盒子
    cubes   =   level.cubes
    # 叶层聚合项
    aggS    =   level.aggS
    aggS   .=   0
    # 多极子数
    nPoles  =   size(aggS, 1)
    # 将本函数内的BLAS设为单线程
    nthds = nthreads()
    BLAS.set_num_threads(1)
    # 避免多线程内存分配问题
    # IanTemps = zeros(CT, nthds, 2)
    # 对盒子循环计算
    @threads for iCube in eachindex(cubes)
        # 盒子信息``
        cube    =   cubes[iCube]
        # 基函数区间
        bfInterval  =   cube.bfInterval
        # 盒子的聚合项
        aggSCube    =   view(aggS, :, :, iCube)
        # 往盒子中心聚合
        @inbounds for n in bfInterval
            In = IVec[n]
            for dr in 1:2, idx in 1:nPoles
                aggSCube[idx, dr]   +=  In*aggSBF[idx, dr, n]
            end
            # @views aggS[:, :, iCube]   .+=   IVec[n] .* aggSBF[:, :, n]
        end #n
    end #iCube
    # 恢复BLAS默认线程以防影响其他多线程函数
    BLAS.set_num_threads(nthds)
    return nothing
end


"""
从子层聚合到本层
tLevel :: 本层
kLevel :: 子层
"""
function agg2HighLevel!(tLevel, kLevel)
    # 本层信息
    cubes   =   tLevel.cubes
    tAggS   =   tLevel.aggS
    # 置零避免累加错误
    tAggS  .=   0
    CT      =   eltype(tAggS)
    # 子层信息
    kAggS   =   kLevel.aggS
    # 从子层到本层的相移
    phaseShiftFromKids  =   tLevel.phaseShiftFromKids
    # 子层到本层的稀疏插值矩阵
    θCSC    =   kLevel.interpWθϕ.θCSC
    ϕCSC    =   kLevel.interpWθϕ.ϕCSC
    # 预分配插值结果内存
    # aggSInterped =  zeros(eltype(tAggS), size(tAggS, 1), size(tAggS, 2))
    
    nthds = nthreads()
    # 对盒子循环
    if nthds > length(cubes) # 线程数大于盒子数，用 BLAS 多线程
        # 分配内存
        aggSInterpedϕ1  =   zeros(CT, size(ϕCSC, 1), size(kAggS, 2))
        aggSInterped1   =   zeros(CT, size(θCSC, 1), size(kAggS, 2))
        for iCube in eachindex(cubes)
            cube    =   cubes[iCube]
            # 对子盒子循环
            # 子盒子数
            nkCube  =   length(cube.kidsInterval)                
            @inbounds for jkCube in 1:nkCube
                # 子盒子id
                kCubeID =   cube.kidsInterval[jkCube]
                # 子盒子在8个子盒子中的id
                kIn8    =   cube.kidsIn8[jkCube]
                
                # 将子盒子的辐射积分进行插值，必须按先 ϕ 后 θ 的顺序
                @views mul!(aggSInterpedϕ1, ϕCSC, kAggS[:, :, kCubeID])
                mul!(aggSInterped1, θCSC, aggSInterpedϕ1)

                # 插值完毕进行子层到本层的相移，累加进父盒子聚合项
                @views tAggS[:,:,iCube]   .+=   phaseShiftFromKids[:,kIn8] .* aggSInterped1
            end # jkCube
        end #iCube
    else # 线程数小于盒子数，用 @threads 多线程
        # 将本函数内的BLAS设为单线程
        BLAS.set_num_threads(1)
        # 分线程分配内存
        aggSInterpedϕs  =   zeros(CT, size(ϕCSC, 1), size(kAggS, 2), nthds)
        aggSInterpeds   =   zeros(CT, size(θCSC, 1), size(kAggS, 2), nthds)
        # 各个线程预分配内存
        @threads for iCube in eachindex(cubes)
            cube    =   cubes[iCube]
            # 线程 id
            tid     =   Threads.threadid()
            aggSInterpedϕ2   =   view(aggSInterpedϕs, :, :, tid)
            aggSInterped2    =   view(aggSInterpeds, :, :, tid)
            # 对子盒子循环
            # 子盒子数
            nkCube  =   length(cube.kidsInterval)
            @inbounds for jkCube in 1:nkCube
                # 子盒子id
                kCubeID =   cube.kidsInterval[jkCube]
                # 子盒子在8个子盒子中的id
                kIn8    =   cube.kidsIn8[jkCube]
                
                # 将子盒子的辐射积分进行插值，必须按先 ϕ 后 θ 的顺序
                # @views aggSInterpedϕ .= ϕCSC*kAggS[:, :, kCubeID]
                # @views aggSInterped  .= θCSC*aggSInterpedϕ
                @views mul!(aggSInterpedϕ2, ϕCSC, kAggS[:, :, kCubeID])
                mul!(aggSInterped2, θCSC, aggSInterpedϕ2)

                # 插值完毕进行子层到本层的相移，累加进父盒子聚合项
                @views tAggS[:,:,iCube]   .+=   phaseShiftFromKids[:,kIn8] .* aggSInterped2

            end # jkCube
        end #iCube
        # 恢复BLAS默认线程以防影响其他多线程函数
        BLAS.set_num_threads(nthds)
    end # if

end # function


"""
从叶层聚合到第 '2' 层
"""
function agg2Level2!(levels, nLevels)
    # 对层循环进行计算
    @inbounds for iLevel in (nLevels - 1):-1:2
        # 本层
        tLevel  =   levels[iLevel]
        # 子层
        kLevel  =   levels[iLevel + 1]
        # 计算
        agg2HighLevel!(tLevel, kLevel)
    end #iLevel
end # function

"""
层内转移
"""
function transOnLevel!(level)
    # 层信息
    cubes   =   level.cubes
    aggS    =   level.aggS
    disaggG =   level.disaggG
    # 置零避免累加错误
    disaggG .=   0
    # 本层的316个转移因子和其索引 OffsetArray 矩阵
    αTrans  =   level.αTrans
    αTransIndex =   level.αTransIndex
    # 对盒子循环
    nthds = nthreads()
    # 对盒子循环
    if nthds > length(cubes) # 线程数大于盒子数，用 BLAS 多线程
        @inbounds for iCube in eachindex(cubes)
            cube    =   cubes[iCube]
            farNeighborIDs = cube.farneighbors
            # 对远亲循环
            for iFarNei in farNeighborIDs
                # 远亲盒子3DID
                farNeiCube3D    =   cubes[iFarNei].ID3D
                # 本盒子相对远亲盒子id
                relative3DID    =   cube.ID3D .- farNeiCube3D
                # 相对 id
                i1d =   αTransIndex[relative3DID[1], relative3DID[2], relative3DID[3]]
                # 转移
                @views disaggG[:, :, iCube]    .+=   αTrans[:, i1d] .* aggS[:, :, iFarNei]
            end # iFarNei
        end #iCube
    else  # 线程数小于盒子数，用 @threads 多线程
        # 将本函数内的BLAS设为单线程，并行采用 @threads 以提高速度@threads
        BLAS.set_num_threads(1)
        @threads for iCube in eachindex(cubes)
            cube    =   cubes[iCube]
            farNeighborIDs = cube.farneighbors
            # 对远亲循环
            @inbounds for iFarNei in farNeighborIDs
                # 远亲盒子3DID
                farNeiCube3D    =   cubes[iFarNei].ID3D
                # 本盒子相对远亲盒子id
                relative3DID    =   cube.ID3D - farNeiCube3D
                # 相对 id
                i1d =   αTransIndex[relative3DID[1], relative3DID[2], relative3DID[3]]
                # 转移
                @views disaggG[:, :, iCube]    .+=   αTrans[:, i1d] .* aggS[:, :, iFarNei]
            end # iFarNei
        end #iCube
        # 恢复BLAS默认线程以防影响其他多线程函数
        BLAS.set_num_threads(nthds)
    end # if
end #function

"""
各层内转移
"""
function transOnLevels!(levels, nLevels)
    @inbounds for iLevel in 2:nLevels
        # 层信息
        level   =   levels[iLevel]
        # 计算
        transOnLevel!(level)
    end #iLevel
end #function


"""
向低层解聚
tLevel :: 本层
kLevel :: 子层
"""
function disagg2KidLevel!(tLevel, kLevel)

    # 本层信息
    cubes   =   tLevel.cubes
    tDisaggG    =   tLevel.disaggG
    # 子层信息
    kDisAggG    =   kLevel.disaggG
    CT          =   eltype(kDisAggG)
    # 从本层到子层的相移
    phaseShift2Kids  =   tLevel.phaseShift2Kids
    # 本层到子层的稀疏反插值矩阵
    θCSCT    =   kLevel.interpWθϕ.θCSCT
    ϕCSCT    =   kLevel.interpWθϕ.ϕCSCT
    # 对盒子循环
    nthds = nthreads()
    # 对盒子循环
    if nthds > length(cubes) # 线程数大于盒子数，用 BLAS 多线程
        # 对盒子循环
        # 分配内存
        disGshifted1    =   zeros(CT, size(tDisaggG, 1), size(tDisaggG, 2))
        disGInterpedθ1  =   zeros(CT, size(θCSCT, 1), size(tDisaggG, 2))
        disGInterped1   =   zeros(CT, size(ϕCSCT, 1), size(tDisaggG, 2))
        @inbounds for iCube in eachindex(cubes)
            cube    =   cubes[iCube]
            @views tCubeDisaggG    =   tDisaggG[:,:,iCube]
            # 对子盒子循环
            # 子盒子数
            nkCube  =   length(cube.kidsInterval)              
            for jkCube in 1:nkCube
                # 子盒子id
                kCubeID =   cube.kidsInterval[jkCube]
                # 子盒子在8个子盒子中的id
                kIn8    =   cube.kidsIn8[jkCube]

                # 进行本层到子层的相移，累加进父盒子聚合项
                @views disGshifted1 .= phaseShift2Kids[:,kIn8] .* tCubeDisaggG
                mul!(disGInterpedθ1, θCSCT, disGshifted1)
                mul!(disGInterped1, ϕCSCT, disGInterpedθ1)
                # @views kdisAggGAnterped = ϕCSCT*(θCSCT*(phaseShift2Kids[:,kIn8] .* tCubeDisaggG))
                kDisAggG[:,:,kCubeID]   .+=  disGInterped1
            end # jkCube
        end #iCube
    else  # 线程数小于盒子数，用 @threads 多线程
        # 将本函数内的BLAS设为单线程
        BLAS.set_num_threads(1)
        # 分线程分配内存
        disGshifteds    =   zeros(CT, size(tDisaggG, 1), size(tDisaggG, 2), nthds)
        disGInterpedθs  =   zeros(CT, size(θCSCT, 1), size(tDisaggG, 2), nthds)
        disGInterpeds   =   zeros(CT, size(ϕCSCT, 1), size(tDisaggG, 2), nthds)
        # 对盒子循环
        @threads for iCube in eachindex(cubes)
            cube    =   cubes[iCube]
            @views tCubeDisaggG    =   tDisaggG[:,:,iCube]
            # 线程 id
            tid     =   Threads.threadid()
            disGshifted2    =   view(disGshifteds, :, :, tid)
            disGInterpedθ2  =   view(disGInterpedθs, :, :, tid)
            disGInterped2   =   view(disGInterpeds, :, :, tid)
            # 对子盒子循环
            # 子盒子数
            nkCube  =   length(cube.kidsInterval)              
            @inbounds for jkCube in 1:nkCube
                # 子盒子id
                kCubeID =   cube.kidsInterval[jkCube]
                # 子盒子在8个子盒子中的id
                kIn8    =   cube.kidsIn8[jkCube]

                # 进行本层到子层的相移，累加进父盒子聚合项
                @views disGshifted2 .= phaseShift2Kids[:,kIn8] .* tCubeDisaggG
                mul!(disGInterpedθ2, θCSCT, disGshifted2)
                mul!(disGInterped2, ϕCSCT, disGInterpedθ2)

                # @views kdisAggGAnterped     =   ϕCSCT*(θCSCT*(phaseShift2Kids[:,kIn8] .* tCubeDisaggG))
                kDisAggG[:,:,kCubeID]     .+=   disGInterped2
            end # jkCube
        end #iCube
        # 恢复BLAS默认线程以防影响其他多线程函数
        BLAS.set_num_threads(nthds)
    end # if

end #function

"""
解聚到叶层
"""
function disagg2LeafLevel!(levels, nLevels)

    # 对层循环进行计算
    @inbounds for iLevel in 2:(nLevels - 1)
        # 本层
        tLevel  =   levels[iLevel]
        # 子层
        kLevel  =   levels[iLevel + 1]
        # 计算
        disagg2KidLevel!(tLevel, kLevel)
    end #iLevel

end #function

"""
在叶层往测试基函数解聚
"""
function disaggOnBF!(level, disaggSBF, ZI)
    CT = eltype(disaggSBF)
    
    # 叶层盒子
    cubes   =   level.cubes
    # 叶层解聚项
    disaggG =   level.disaggG
    # 多极子数
    nPoles  =   size(disaggG, 1)
    # 常量
    JK_0η::CT   =   Params.JK_0*η_0
    # 将本函数内的BLAS设为单线程
    nthds = nthreads()
    BLAS.set_num_threads(1)
    # 避免多线程内存分配问题
    ZInTemps = zeros(CT, nthds)
    # 对盒子循环计算
    @threads for iCube in eachindex(cubes)
        # 盒子信息
        cube    =   cubes[iCube]
        tid     =   Threads.threadid()
        # 基函数区间
        bfInterval  =   cube.bfInterval
        # 该盒子解聚项
        @views disaggGCube =  disaggG[:,:, iCube]
        # 往基函数解聚
        @inbounds for n in bfInterval
            ZInTemps[tid]   = 0
            ZInTemp     =   ZInTemps[tid]
            for idx in 1:nPoles
                ZInTemp += disaggSBF[idx, 1, n] * disaggGCube[idx, 1] + disaggSBF[idx, 2, n] * disaggGCube[idx, 2]
            end
            ZInTemp *=  JK_0η
            # 计算完毕写入数据
            ZI[n]   +=  ZInTemp
        end #n
    end #iCube
    # 恢复BLAS默认线程以防影响其他多线程函数
    BLAS.set_num_threads(nthds)
    return nothing
end

"""
计算远区矩阵向量乘积
"""
function calZfarI!(Zopt::MLMFAIterator{ZT, MT}, IVec::AbstractArray{T}; setzero = true) where {ZT, T<:Number, MT<:Vector}

    # 计算前置零
    setzero && fill!(Zopt.ZI, zero(T))

    # 基函数聚合到叶层
    aggOnBF!(Zopt.leafLevel, Zopt.aggSBF,  IVec)
    # 聚合到2层
    agg2Level2!(Zopt.levels, Zopt.nLevels)
    # 层间转移
    transOnLevels!(Zopt.levels, Zopt.nLevels)
    # 解聚到叶层
    disagg2LeafLevel!(Zopt.levels, Zopt.nLevels)
    # 解聚到基函数
    disaggOnBF!(Zopt.leafLevel, Zopt.disaggSBF, Zopt.ZI)
    
    return Zopt.ZI
end
