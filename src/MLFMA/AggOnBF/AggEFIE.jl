
"""
获取几何信息数组的区间，针对普通 Vector 和 OffsetVector 分别派发
"""
function getGeosInterval(geosInfo::T) where {T<:AbstractVector}
    nhexa   =   length(geosInfo)
    return 1:nhexa
end

"""
获取几何信息数组的区间，针对普通 Vector 和 OffsetVector 分别派发
"""
function getGeosInterval(geosInfo::T) where {T<:OffsetVector}
    nhexa   =   length(geosInfo)
    st      =   (eachindex(geosInfo).offset) + 1
    st:(st - 1 + nhexa)
end

"""
计算某层聚合项, 输入为三角形信息和 RWG 基函数信息
"""
function aggSBFOnLevelEFIE(level, trianglesInfo::Vector{TriangleInfo{IT, FT}}, 
    rwgsInfo) where {IT<:Integer, FT<:Real}
    # 层采样点
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # poles索引
    polesIndices    =   eachindex(polesr̂sθsϕs)
    # 采样多极子数量
    nPoles  =   polesIndices.stop
    # 预分配内存
    nrwg    =   length(rwgsInfo)
    aggSBF  =   zeros(Complex{FT}, nPoles, 2, nrwg)
    disaggSBF   =   similar(aggSBF)

    aggSBFOnLevelEFIE!(aggSBF, disaggSBF, level, trianglesInfo, eltype(rwgsInfo))

    return aggSBF, disaggSBF
end

"""
计算某层聚合项, 输入为三角形信息和 RWG 基函数信息
"""
function aggSBFOnLevelEFIE!(aggSBF, disaggSBF, level, trianglesInfo::Vector{TriangleInfo{IT, FT}}, 
    ::Type{BFT}) where {IT<:Integer, FT<:Real, BFT<: RWG}
    CT  =   Complex{FT}
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
    ntri = length(trianglesInfo)

    # Progress Meter
    pmeter  =   Progress(nCubes, "Aggregating on RWG (EFIE)...")
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
                    aggSθ   =   zero(CT)
                    aggSϕ   =   zero(CT)
                    # 对高斯求积点循环
                    for gi in 1:GQPNTri
                        # 公用的 指数项和权重边长
                        expWlntemp =   exp(JK_0*(poler̂θϕ.r̂ ⋅ cubeC2rgs[:,gi]))*(weightTridiv2[gi]*ln)
                        # 在 θϕ 方向累加
                        @views aggSθ += (poler̂θϕ.θhat ⋅ ρs[:, gi])*expWlntemp
                        @views aggSϕ += (poler̂θϕ.ϕhat ⋅ ρs[:, gi])*expWlntemp
                    end # gi 
                    # 将结果写入目标数组
                    aggSBF[iPole, 1, n] +=  aggSθ
                    aggSBF[iPole, 2, n] +=  aggSϕ
                end # iPole
            end # ni
        end #iTri
        # 更新进度条
        next!(pmeter)
    end #iCube

    # 解聚因子计算
    disaggSBF  .=   conj(aggSBF)

    return nothing
end


"""
计算某层聚合项, 输入为四面体信息和 SWG 基函数信息
"""
function aggSBFOnLevel(level::LT, geosInfo::AbstractVector{VT}, 
    bfsInfo::AbstractVector{BFT}) where {LT<:LevelInfo, VT<:VolumeCellType, BFT<:BasisFunctionType}
    # 层采样点
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # poles索引
    polesIndices    =   eachindex(polesr̂sθsϕs)
    # 采样多极子数量
    nPoles  =   polesIndices.stop
    # 预分配内存
    nbf     =   length(bfsInfo)
    FT      =   typeof(level.cubeEdgel)
    aggSBF  =   zeros(Complex{FT}, nPoles, 2, nbf)
    disaggSBF   =   zeros(Complex{FT}, nPoles, 2, nbf)
    # 计算
    aggSBFOnLevel!(aggSBF, disaggSBF, level, geosInfo, eltype(bfsInfo))

    return aggSBF, disaggSBF
end



"""
计算某层聚合项, 输入为四面体信息和 SWG 基函数信息
"""
function aggSBFOnLevel!(aggSBF, disaggSBF, level, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, 
    ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:SWG}
    # 层采样点
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # poles索引
    polesIndices    =   eachindex(polesr̂sθsϕs)
    # 本层盒子信息
    cubes   =   level.cubes
    nCubes  =   eachindex(cubes).stop
    # 四面体高斯求积权重
    weightTetradiv3   =   TetraGQInfo.weight / 3
    # 常数
    JK_0 = Params.JK_0
    # 判断体电流的离散方式
    discreteJ::Bool = SimulationParams.discreteVar == "J"
    # 是否为偏置数组
    geoInterval =   getGeosInterval(tetrasInfo)
    # Progress Meter
    pmeter  =   Progress(nCubes, "Aggregating on SWG (EFIE)...")
    # 对盒子循环计算
    @threads for iCube in eachindex(cubes)
        
        # 盒子
        cube    =   cubes[iCube]
        # 盒子中心
        cubeCenter  =   cube.center

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的四面体id
        cubeTetrasID    =   cube.geoIDs
        # for (i, swg) in enumerate(view(swgsInfo, cubeBFinterval))
        #     cubeTetrasID[(2i-1):(2i)]  .=  swg.inGeo
        # end
        # # 排序并剔除冗余元素
        # unique!(sort!(cubeTetrasID))
        # 防止边界元基函数的负部，编号为 0 的盒子出现在索引中
        cubeTetrasID[1] == 0 && popfirst!(cubeTetrasID)
        # 对盒子内四面体循环
        for iTetra in 1:length(cubeTetrasID)
            it = cubeTetrasID[iTetra]
            # 混合网格不是四面体的id则跳过
            !(it in geoInterval) && continue
            # 四面体
            tetra   =   tetrasInfo[it]
            # 高斯求积点
            rgs     =   getGQPTetra(tetra)
            # 盒子中心到求积点向量
            cubeC2rgs   =   zero(rgs)
            for gi in 1:GQPNTetra
                cubeC2rgs[:, gi]   .=   view(rgs, :, gi) .- cubeCenter
            end
            # 介质对比度
            κₜ  =   tetra.κ

            # 对四面体上的基函数循环
            for ni in 1:4
                # 基函数编号
                n   =   tetra.inBfsID[ni]
                # 基函数不在该盒子的基函数区间则跳过
                !(n in cubeBFinterval) && continue
                # arean
                arean   =   tetra.facesArea[ni]
                # ρs
                ρs      =   zero(rgs)
                for gi in 1:GQPNTetra
                    ρs[:, gi]   .=   view(rgs, :, gi) .- view(tetra.vertices, :, ni)
                end

                # 对多极子循环计算
                for iPole in polesIndices
                    # 该多极子
                    poler̂θϕ =   polesr̂sθsϕs[iPole]
                    # 聚合项初始化
                    aggSθ   =   zero(CT)
                    aggSϕ   =   zero(CT)
                    # 对高斯求积点循环
                    for gi in 1:GQPNTetra
                        # 公用的 指数项和权重边长
                        expWareantemp =   exp(JK_0*(poler̂θϕ.r̂ ⋅ view(cubeC2rgs, :, gi)))*(weightTetradiv3[gi]*arean)
                        # 在 θϕ 方向累加
                        @views aggSθ += (poler̂θϕ.θhat ⋅ ρs[:, gi])*expWareantemp
                        @views aggSϕ += (poler̂θϕ.ϕhat ⋅ ρs[:, gi])*expWareantemp
                    end # gi 
                    # 将结果写入目标数组
                    if  discreteJ
                        aggSBF[iPole, 1, n]     += aggSθ
                        aggSBF[iPole, 2, n]     += aggSϕ
                    else
                        aggSBF[iPole, 1, n]     += κₜ*aggSθ
                        aggSBF[iPole, 2, n]     += κₜ*aggSϕ
                    end
                    disaggSBF[iPole, 1, n]  += conj(aggSθ)
                    disaggSBF[iPole, 2, n]  += conj(aggSϕ)
                end # iPole
            end # ni
        end #iTetra
        # 更新进度条
        next!(pmeter)
    end #iCube

    return nothing
end


"""
计算某层聚合项, 输入为四面体信息和 PWC 基函数信息
"""
function aggSBFOnLevel!(aggSBF, disaggSBF, level, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, 
    ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:PWC}
    # 层采样点
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # poles索引
    polesIndices    =   eachindex(polesr̂sθsϕs)
    # 本层盒子信息
    cubes   =   level.cubes
    nCubes  =   eachindex(cubes).stop
    # 网格区间信息
    ntetra  =   length(tetrasInfo)
    isoffset    =   isa(tetrasInfo, OffsetVector)
    geoInterval =   begin
        isoffset ? begin
            st  =   (eachindex(tetrasInfo).offset) + 1
            st:(st - 1 + ntetra)
        end : begin
            1:ntetra
        end
    end
    # 四面体高斯求积权重
    weightTetra   =   TetraGQInfo.weight
    # 常数
    JK_0 = Params.JK_0
    # 判断体电流的离散方式
    discreteJ::Bool = SimulationParams.discreteVar === "J"
    # Progress Meter
    pmeter  =   Progress(nCubes, "Aggregating on PWC (EFIE)...")
    # 对盒子循环计算
    @threads for iCube in eachindex(cubes)
        
        # 盒子
        cube    =   cubes[iCube]
        # 盒子中心
        cubeCenter  =   cube.center

        # 找出对应的四面体id
        cubeTetrasID    =   cube.geoIDs
        # 防止边界元基函数的负部，编号为 0 的盒子出现在索引中
        cubeTetrasID[1] == 0 && popfirst!(cubeTetrasID)
        # 对盒子内四面体循环
        for iTetra in 1:length(cubeTetrasID)
            it = cubeTetrasID[iTetra]
            # 混合网格不是四面体的id则跳过
            !(it in geoInterval) && continue
            # 四面体
            tetra   =   tetrasInfo[it]
            # 高斯求积点
            rgs     =   getGQPTetra(tetra)
            # 盒子中心到求积点向量
            cubeC2rgs   =   zero(rgs)
            for gi in 1:GQPNTetra
                cubeC2rgs[:, gi]   .=   view(rgs, :, gi) .- cubeCenter
            end
            # 介质对比度
            κₜ  =   tetra.κ
            # 体积
            Vₜ  =   tetra.volume

            # 对四面体上的基函数循环
            for ni in 1:3
                # 基函数编号
                n   =   tetra.inBfsID[ni]

                # 对多极子循环计算
                for iPole in polesIndices
                    # 该多极子
                    poler̂θϕ =   polesr̂sθsϕs[iPole]
                    # 聚合项初始化
                    aggSθ   =   zero(CT)
                    aggSϕ   =   zero(CT)
                    # 对高斯求积点循环
                    for gi in 1:GQPNTetra
                        # 公用的 指数项和权重边长
                        expWVtemp =   exp(JK_0*(poler̂θϕ.r̂ ⋅ view(cubeC2rgs, :, gi)))*(weightTetra[gi]*Vₜ)
                        # 在 θϕ 方向累加
                        @views aggSθ += (poler̂θϕ.θhat[ni])*expWVtemp
                        @views aggSϕ += (poler̂θϕ.ϕhat[ni])*expWVtemp
                    end # gi 
                    # 将结果写入目标数组
                    if  discreteJ
                        aggSBF[iPole, 1, n]     += aggSθ
                        aggSBF[iPole, 2, n]     += aggSϕ
                    else
                        aggSBF[iPole, 1, n]     += κₜ*aggSθ
                        aggSBF[iPole, 2, n]     += κₜ*aggSϕ
                    end
                    disaggSBF[iPole, 1, n]  += conj(aggSθ)
                    disaggSBF[iPole, 2, n]  += conj(aggSϕ)
                end # iPole
            end # ni
        end #iTetra
        # 更新进度条
        next!(pmeter)
    end #iCube

    return nothing
end



"""
计算某层聚合项, 输入为六面体信息和 PWC 基函数信息
"""
function aggSBFOnLevel!(aggSBF, disaggSBF, level, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}, 
    ::Type{BFT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:PWC}
    # 层采样点
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # poles索引
    polesIndices    =   eachindex(polesr̂sθsϕs)

    # 本层盒子信息
    cubes   =   level.cubes
    nCubes  =   eachindex(cubes).stop
    # 四面体高斯求积权重
    weightHexa   =   HexaGQInfo.weight
    # 判断体电流的离散方式
    discreteJ::Bool = SimulationParams.discreteVar == "J"
    # 几何信息索引区间
    geoInterval =   getGeosInterval(hexasInfo)
    # 常数
    JK_0 = Params.JK_0
    # Progress Meter
    pmeter  =   Progress(nCubes, "Aggregating on PWC (EFIE)...")
    # 对盒子循环计算
    @threads for iCube in eachindex(cubes)
        
        # 盒子
        cube    =   cubes[iCube]
        # 盒子中心
        cubeCenter  =   cube.center

        # 找出对应的四面体id
        cubeHexasID    =   cube.geoIDs
        # 防止边界元基函数的负部，编号为 0 的盒子出现在索引中
        cubeHexasID[1] == 0 && popfirst!(cubeHexasID)
        # 对盒子内六面体循环
        for iHexa in 1:length(cubeHexasID)
            it  =   cubeHexasID[iHexa]
            # 混合网格时跳过非本类型网格的id
            !(it in geoInterval) && continue
            # 六面体
            hexa    =   hexasInfo[it]
            # 高斯求积点
            rgs     =   getGQPHexa(hexa)
            # 盒子中心到求积点向量
            cubeC2rgs   =   zero(rgs)
            for gi in 1:GQPNHexa
                cubeC2rgs[:, gi]   .=   view(rgs, :, gi) .- cubeCenter
            end
            # 介质对比度
            κₜ  =   hexa.κ
            # 体积
            Vₜ  =   hexa.volume

            # 对六面体上的基函数循环
            for ni in 1:3
                # 基函数编号
                n   =   hexa.inBfsID[ni]

                # 对多极子循环计算
                for iPole in polesIndices
                    # 该多极子
                    poler̂θϕ =   polesr̂sθsϕs[iPole]
                    # 聚合项初始化
                    aggSθ   =   zero(CT)
                    aggSϕ   =   zero(CT)
                    # 对高斯求积点循环
                    for gi in 1:GQPNHexa
                        # 公用的 指数项和权重边长
                        expWVtemp =   exp(JK_0*(poler̂θϕ.r̂ ⋅ view(cubeC2rgs, :, gi)))*(weightHexa[gi]*Vₜ)
                        # 在 θϕ 方向累加
                        @views aggSθ += (poler̂θϕ.θhat[ni])*expWVtemp
                        @views aggSϕ += (poler̂θϕ.ϕhat[ni])*expWVtemp
                    end # gi 
                    # 将结果写入目标数组
                    if  discreteJ
                        aggSBF[iPole, 1, n]     += aggSθ
                        aggSBF[iPole, 2, n]     += aggSϕ
                    else
                        aggSBF[iPole, 1, n]     += κₜ*aggSθ
                        aggSBF[iPole, 2, n]     += κₜ*aggSϕ
                    end
                    disaggSBF[iPole, 1, n]  += conj(aggSθ)
                    disaggSBF[iPole, 2, n]  += conj(aggSϕ)
                end # iPole
            end # ni
        end #iHexa
        # 更新进度条
        next!(pmeter)
    end #iCube

    return nothing
end

"""
计算某层聚合项, 输入为四面体信息和 SWG 基函数信息
"""
function aggSBFOnLevel!(aggSBF, disaggSBF, level, hexasInfo::AbstractVector{VT}, 
    ::Type{BFT}) where {VT<:HexahedraInfo, BFT<:RBF}
    CT  =   Precision.CT
    # 层采样点
    polesr̂sθsϕs =   level.poles.r̂sθsϕs
    # poles索引
    polesIndices    =   eachindex(polesr̂sθsϕs)
    # 本层盒子信息
    cubes   =   level.cubes
    nCubes  =   eachindex(cubes).stop
    # 四面体高斯求积权重
    weightHexa  =   HexaGQInfo.weight
    # 常数
    JK_0 = Params.JK_0
    # 判断体电流的离散方式
    discreteJ::Bool = SimulationParams.discreteVar == "J"
    # 几何信息索引区间
    geoInterval =   getGeosInterval(hexasInfo)
    
    # Progress Meter
    pmeter  =   Progress(nCubes, "Aggregating on RBF (EFIE)...")
    # 对盒子循环计算
    @threads for iCube in eachindex(cubes)
        
        # 盒子
        cube    =   cubes[iCube]
        # 盒子中心
        cubeCenter  =   cube.center

        # 盒子里的基函数区间
        cubeBFinterval  =   cube.bfInterval
        # 找出对应的四面体id
        cubeHexasID    =   cube.geoIDs
        # for (i, swg) in enumerate(view(rbfsInfo, cubeBFinterval))
        #     cubeHexasID[(2i-1):(2i)]  .=  swg.inGeo
        # end
        # # 排序并剔除冗余元素
        # unique!(sort!(cubeHexasID))
        # 防止边界元基函数的负部，编号为 0 的盒子出现在索引中
        cubeHexasID[1] == 0 && popfirst!(cubeHexasID)
        # 对盒子内四面体循环
        for iHexa in 1:length(cubeHexasID)
            it = cubeHexasID[iHexa]
            # 跳过混合网格时非本类型网格的 id
            !(it in geoInterval) && continue
            # 四面体
            hexa   =   hexasInfo[it]
            # 高斯求积点
            rgs     =   getGQPHexa(hexa)
            # 盒子中心到求积点向量
            cubeC2rgs   =   zero(rgs)
            for gi in 1:GQPNHexa
                cubeC2rgs[:, gi]   .=   view(rgs, :, gi) .- cubeCenter
            end
            # 介质对比度
            κₜ  =   hexa.κ

            # 对六面体上的基函数循环
            for ni in 1:6
                freeVns     =   getFreeVns(hexa, ni)
                # 基函数编号
                n   =   hexa.inBfsID[ni]
                # 基函数不在该盒子的基函数区间则跳过
                !(n in cubeBFinterval) && continue
                # arean
                arean   =   hexa.facesArea[ni]
                # ρs
                ρs      =   zero(rgs)
                for gi in 1:GQPNHexa
                    # 计算场六面体的 “自由端” id
                    idn3D   =   GQ1DID2GQ3DIDVector[gi]
                    # 计算源六面体的 “自由端” id
                    idn     =   getFreeVIDFromGQ3DID(idn3D, ni)
                    # ρs
                    ρs[:, gi]   .=   view(rgs, :, gi) .- view(freeVns, :, idn)
                end

                # 对多极子循环计算
                for iPole in polesIndices
                    # 该多极子
                    poler̂θϕ =   polesr̂sθsϕs[iPole]
                    # 聚合项初始化
                    aggSθ   =   zero(CT)
                    aggSϕ   =   zero(CT)
                    # 对高斯求积点循环
                    for gi in 1:GQPNHexa
                        # 公用的 指数项和权重边长
                        expWareantemp =   exp(JK_0*(poler̂θϕ.r̂ ⋅ view(cubeC2rgs, :, gi)))*(weightHexa[gi]*arean)
                        # 在 θϕ 方向累加
                        @views aggSθ += (poler̂θϕ.θhat ⋅ ρs[:, gi])*expWareantemp
                        @views aggSϕ += (poler̂θϕ.ϕhat ⋅ ρs[:, gi])*expWareantemp
                    end # gi 
                    # 将结果写入目标数组
                    if  discreteJ
                        aggSBF[iPole, 1, n]     += aggSθ
                        aggSBF[iPole, 2, n]     += aggSϕ
                    else
                        aggSBF[iPole, 1, n]     += κₜ*aggSθ
                        aggSBF[iPole, 2, n]     += κₜ*aggSϕ
                    end
                    disaggSBF[iPole, 1, n]  += conj(aggSθ)
                    disaggSBF[iPole, 2, n]  += conj(aggSϕ)
                end # iPole
            end # ni
        end #iHexa
        # 更新进度条
        next!(pmeter)
    end #iCube
    println("ok!")

    return nothing
end