"""
    aggSBFOnLevelMFIE(level, trianglesInfo::Vector{TriangleInfo{IT, FT}}, 
    bfsInfo::Vector{RWG{IT, FT}}) where {IT<:Integer, FT<:Real}

计算某层采用 MFIE 时在三角形上的 RWG 基函数的辐射函数 `aggSBF` 、配置函数 `disaggSBF`。
"""
function aggSBFOnLevelMFIE(level, trianglesInfo::Vector{TriangleInfo{IT, FT}}, 
    bfsInfo::Vector{RWG{IT, FT}}) where {IT<:Integer, FT<:Real}
    CT  =   Complex{FT}
    # 层采样点
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # poles索引
    polesIndices    =   eachindex(polesr̂sθsϕs)
    # 采样多极子数量
    nPoles  =   polesIndices.stop
    # 预分配内存
    nbf    =   length(bfsInfo)
    aggSBF      =   zeros(CT, nPoles, 2, nbf)
    disaggSBF   =   zeros(CT, nPoles, 2, nbf)
    # 计算
    aggSBFOnLevelMFIE!(aggSBF, disaggSBF, level, trianglesInfo, eltype(bfsInfo))

    return aggSBF, disaggSBF
end

using MoM_Basics:acrossbcrossc
"""
    aggSBFOnLevelMFIE!(aggSBF, disaggSBF, level, trianglesInfo::Vector{TriangleInfo{IT, FT}}, 
    ::Type{BFT}; setzero = true) where {IT<:Integer, FT<:Real, BFT<:RWG}

计算某层采用 MFIE 时在三角形上的 RWG 基函数的辐射函数 `aggSBF` 、配置函数 `disaggSBF`。
"""
function aggSBFOnLevelMFIE!(aggSBF, disaggSBF, level, trianglesInfo::Vector{TriangleInfo{IT, FT}}, 
    ::Type{BFT}; setzero = true) where {IT<:Integer, FT<:Real, BFT<:RWG}
    CT  =   Complex{FT}
    setzero && begin 
        fill!(aggSBF, 0)
        fill!(disaggSBF, 0)
    end
    # 层采样点
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # poles索引
    polesIndices    =   eachindex(polesr̂sθsϕs)
    # 本层盒子信息
    cubes   =   level.cubes
    nCubes  =   eachindex(cubes).stop
    # 三角形高斯求积权重
    weightTridiv2   =   TriGQInfo.weight / 2
    # 常数
    JK_0 = Params.JK_0
    CT0  = zero(CT)
    ntri = length(trianglesInfo)

    # Progress Meter
    pmeter  =   Progress(nCubes, "Aggregating on RWG (MFIE)...")
    # 对盒子循环计算
    @threads for iCube in eachindex(cubes)
        # 盒子
        cube    =   cubes[iCube]
        # 盒子中心
        cubeCenter  =   cube.center

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的三角形id
        cubeTriID       =   cube.geoIDs
        # 排序并剔除冗余元素
        # unique!(sort!(cubeTriID))
        # 对盒子内三角形循环
        for iTri in 1:length(cubeTriID)
            it = cubeTriID[iTri]
            # 超出区间跳过
            it > ntri && continue
            # 三角形
            tri =   trianglesInfo[it]
            # 面外法向量
            n̂   =   tri.facen̂
            # 高斯求积点
            rgs =   getGQPTri(tri)
            # 盒子中心到求积点向量
            cubeC2rgs   =   zero(rgs)
            for gi in 1:GQPNTri
                cubeC2rgs[:, gi]   .=   view(rgs, :, gi) .- cubeCenter
            end

            # 对三角形上的基函数循环
            for ni in 1:3
                # 基函数编号
                n   =   tri.inBfsID[ni]
                # 基函数不在该盒子的基函数区间则跳过
                !(n in cubeBFinterval) && continue
                # ln
                ln  =   tri.edgel[ni]
                # ρs
                ρs  =   zero(rgs)
                @views for gi in 1:GQPNTri
                    ρs[:, gi]   .=   rgs[:, gi] .- tri.vertices[:,ni]
                end

                # 对多极子循环计算
                for iPole in polesIndices
                    # 该多极子
                    poler̂θϕ =   polesr̂sθsϕs[iPole]
                    # 聚合项初始化
                    aggSθ   =   CT0
                    aggSϕ   =   CT0
                    disaggSθ    =   CT0
                    disaggSϕ    =   CT0
                    # 对高斯求积点循环
                    for gi in 1:GQPNTri
                        ρj   =   ρs[:, gi]
                        # 公用的 指数项和权重边长
                        expWlntemp  =   exp(JK_0*(poler̂θϕ.r̂ ⋅ cubeC2rgs[:,gi]))*(weightTridiv2[gi]*ln)
                        # 在 θϕ 方向累加
                        @views aggSθ += (poler̂θϕ.θhat ⋅ ρj)*expWlntemp
                        @views aggSϕ += (poler̂θϕ.ϕhat ⋅ ρj)*expWlntemp
                        # ρj × n̂ × k̂
                        ρcn̂ck̂   =   acrossbcrossc(ρj, n̂, poler̂θϕ.r̂)
                        @views disaggSθ += (poler̂θϕ.θhat ⋅ ρcn̂ck̂)*expWlntemp'
                        @views disaggSϕ += (poler̂θϕ.ϕhat ⋅ ρcn̂ck̂)*expWlntemp'
                    end # gi 
                    # 将结果写入目标数组
                    aggSBF[iPole, 1, n]     +=  aggSθ
                    aggSBF[iPole, 2, n]     +=  aggSϕ
                    disaggSBF[iPole, 1, n]  +=  disaggSθ
                    disaggSBF[iPole, 2, n]  +=  disaggSϕ
                end # iPole
            end # ni
        end #iTri
        # 更新进度条
        next!(pmeter)
    end #iCube

    return nothing
end
