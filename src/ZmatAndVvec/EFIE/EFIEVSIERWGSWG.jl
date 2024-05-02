
# 计算时用到的常数

"""
计算三角形和四面体上相关的 12 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较远的情况
输入
tris    ::  TriangleInfo,  源三角形面体
tetrat  ::  TetrahedraInfo, 场四面体
"""
function EFIEOnRWGSWG(trit::TriangleInfo{IT, FT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的小数组
    Zts =   zeros(CT, 3, 4)
    Zst =   zeros(CT, 4, 3)
    
    # 场源求积点
    rgt =   getGQPTri(trit)
    rgs =   getGQPTetra(tetras)
    # 源介质对比度
    κs  =   tetras.κ

    # 常数项
    JKη_0div4π  =   Params.JKη_0/4/π
    divk²       =   1/Params.k²

    # 源四面体包含的面
    facess  =   tetras.faces

    ## 求g乘以权重
    # 对源求积点循环
    # 预分配内存
    gw  =   zero(MMatrix{GQPNTri, GQPNTetra, Complex{FT}})
    @inbounds for gj in 1:GQPNTetra
        # 源高斯求积点
        @views rgj  =   rgs[:, gj]
        # 对场求积点循环
        for gi in 1:GQPNTri
            # 场高斯求积点
            @views rgi  =  rgt[:, gi]
            # g乘以权重
            gw[gi, gj]  =   greenfunc(rgi, rgj)*(TriGQInfo.weight[gi]*TetraGQInfo.weight[gj])
        end # for gi
    end #for gj

    # 预分配数值以加速
    ρnj     =   zero(MVec3D{FT})
    ρmi     =   zero(MVec3D{FT})
    # 对源基函数循环
    @inbounds for ni in 1:4
        # 源基函数带符号面积、自由端
        arean   =   tetras.facesArea[ni]
        freeVn  =   tetras.vertices[:, ni]
        # 第 ni 个面
        faces   =   facess[ni]
        # 介质对比度变化量
        δκn     =   faces.δκ
        # 是否为半基函数？
        isbdn   =   faces.isbd
        # 对场基函数循环
        for mi in 1:3
            # 所在基函数id
            m       =   trit.inBfsID[mi]
            # 判断边是不是基函数（RWG 边缘不算因此要跳过）
            m  == 0 && continue
            # 场基函数序号、n带符号边长、自由端
            lm      =   trit.edgel[mi]
            freeVm  =   trit.vertices[:, mi]
            # 面积乘积
            lman    =   lm*arean

            # 初始化阻抗矩阵元
            Zmn     =   zero(CT)
            Znm     =   zero(CT)

            # 计算 Fsv 项
            Fsv  =   zero(CT)
            # 非奇异直接高斯求积
            for gj in 1:GQPNTetra
                # 源高斯求积点
                @views rgj =   rgs[:, gj]
                # ρn
                ρnj    .=   rgj .- freeVn
                # 对场求积点循环
                for gi in 1:GQPNTri
                    # 场高斯求积点
                    @views rgi =   rgt[:, gi]
                    # ρm
                    ρmi    .=   rgi .- freeVm
                    # 计算结果
                    Fsv  +=  ((ρmi ⋅ ρnj) / 6 - divk²) * gw[gi, gj]                    
                end #for gi
            end #for gj
            Fsv *=  lman
            # 累加 Zmn\Znm
            JKη_0div4πFsv    =   JKη_0div4π*Fsv
            Zmn    +=  κs * JKη_0div4πFsv
            Znm    +=       JKη_0div4πFsv
            ## Fss 项只有满足以下条件时需要计算
            if isbdn || !iszero(δκn)
                # 初始化
                Fss     =   zero(CT)
                # 两个面的距离
                r   =   norm(trit.center .- view(mean(faces.vertices, dims = 2), :, 1))
                # 足够远不处理奇异性
                if r > Params.Rsglr
                    for gj in 1:GQPNTri
                        # 源高斯求积点
                        rgj  =   getGQPTri(faces, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNTri
                            # 场高斯求积点
                            @views rgi =   rgt[:, gi]
                            # g乘以权重
                            Fss  +=   greenfunc(rgi, rgj)*(TriGQInfo.weight[gi]*TriGQInfo.weight[gj])
                        end # for gi
                    end #for gj
                # 过小需要处理奇异性
                else
                    for gj in 1:GQPNTriSglr
                        # 源高斯求积点
                        rgj =   getGQPTriSglr(faces, gj)
                        # IgS  =   ∫ g(R) dS 的奇异性计算
                        Igw =   faceSingularityIg(rgj, trit, trit.area, trit.facen̂)*TriGQInfoSglr.weight[gj]
                        Fss    +=   Igw
                        # 对场求积点循环
                    end #for gj
                    Fss /= trit.area
                end
                # Fss项的系数
                Fss *=  lm*abs(arean)

                # 写入数据
                temp = JKη_0div4π*divk²*Fss
                !iszero(δκn)    &&  begin Zmn += δκn*temp;  end
                isbdn           &&  begin Znm += temp;      end
            end
            # 写入数据
            Zts[mi, ni] =   Zmn
            Zst[ni, mi] =   Znm

        end #for mi
    end # for ni 
    return Zts, Zst
end

"""
计算三角形和四面体上相关的 12 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较远的情况
输入
tris    ::  TriangleInfo,  源三角形面体
tetrat  ::  TetrahedraInfo, 场四面体
"""
function EFIEOnRWGSWG(tetrat::TetrahedraInfo{IT, FT, CT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    Zst, Zts = EFIEOnRWGSWG(tris, tetrat)
    return Zts, Zst
end

"""
计算三角形和四面体上相关的 12 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较近的情况
输入
tris    ::  TriangleInfo,  源三角形面体
tetrat  ::  TetrahedraInfo, 场四面体
"""
function EFIEOnNearRWGSWG(trit::TriangleInfo{IT, FT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts =   zeros(CT, 3, 4)
    Zst =   zeros(CT, 4, 3)
    
    # 场源求积点
    rgt =   getGQPTriSglr(trit)
    rgs =   getGQPTetraSglr(tetras)
    # 源介质对比度
    κs  =   tetras.κ

    # 常数项
    JKη_0div4π  =   Params.JKη_0/4/π
    divk²       =   1/Params.k²

    # 在源积分点计算场三角得到的 Ig, Ivecg
    Igs     =   zero(MVector{GQPNTetraSglr, Complex{FT}})
    Ivecgs  =   zero(MMatrix{3, GQPNTetraSglr, Complex{FT}})
    # 循环计算求积点处的奇异项
    for gj in 1:GQPNTetraSglr
        Ig, Ivecg   =   faceSingularityIgIvecg(rgs[:, gj], trit, abs(trit.area), trit.facen̂)
        Igs[gj]         =   Ig
        Ivecgs[:, gj]  .=   Ivecg
    end
    # 源四面体包含的面
    facess  =   tetras.faces
    # 预分配数值以加速
    ρnjm    =   zero(MVec3D{FT})
    ρnj     =   zero(MVec3D{FT})
    # 对源基函数循环
    @inbounds for ni in 1:4
        # 源基函数带符号面积、自由端
        arean   =   tetras.facesArea[ni]
        freeVn  =   tetras.vertices[:, ni]
        # 第 ni 个面
        faces   =   facess[ni]
        # 介质对比度变化量
        δκn     =   faces.δκ
        # 是否为半基函数？
        isbdn   =   faces.isbd
        # 所在基函数id
        n       =   tetras.inBfsID[ni]
        # 对场基函数循环
        for mi in 1:3
            # 所在基函数id
            m       =   trit.inBfsID[mi]
            # 判断边是不是基函数（RWG 边缘不算因此要跳过）
            m   == 0 && continue
            # 场基函数序号、n带符号边长、自由端
            lm      =   trit.edgel[mi]
            freeVm  =   trit.vertices[:, mi]
            # 面积乘积
            lman    =   lm*arean

            # 初始化阻抗矩阵元
            Zmn     =   zero(CT)
            Znm     =   zero(CT)

            # 计算 Fsv 项
            Fsv  =   zero(CT)
            # 对场求积点循环
            for gj in 1:GQPNTetraSglr
                # 场高斯求积点
                @views rgj =   rgs[:, gj]
                # ρn
                ρnj    .=   rgj .- freeVn
                # ρnim = rgj - freeVm
                ρnjm    =   rgj .- freeVm
                # 计算结果
                Fsv   +=  ((ρnj ⋅ ρnjm / 6 - divk²)*Igs[gj] - ( ρnj ⋅ view(Ivecgs, :, gj) / 6 )) * TetraGQInfoSglr.weight[gj]/trit.area

            end #for gi
            Fsv *=  lman
            # 累加 Zmn\Znm
            JKη_0div4πFsv    =   JKη_0div4π*Fsv
            Zmn    +=  κs * JKη_0div4πFsv
            Znm    +=       JKη_0div4πFsv
            ## Fss 项只有满足以下条件时需要计算
            if isbdn || !iszero(δκn)
                # 初始化
                Fss     =   zero(CT)
                # 两个面的距离
                r   =   norm(trit.center .- view(mean(faces.vertices, dims = 2), :, 1))
                # 足够远不处理奇异性
                if r > Params.Rsglr
                    for gj in 1:GQPNTri
                        # 源高斯求积点
                        rgj  =   getGQPTri(faces, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNTri
                            # 场高斯求积点
                            @views rgi =   rgt[:, gi]
                            # g乘以权重
                            Fss  +=   greenfunc(rgi, rgj)*(TriGQInfo.weight[gi]*TriGQInfo.weight[gj])
                        end # for gi
                    end #for gj
                # 过小需要处理奇异性
                else
                    for gj in 1:GQPNTriSglr
                        # 源高斯求积点
                        rgj =   getGQPTriSglr(faces, gj)
                        # IgS  =   ∫ g(R) dS 的奇异性计算
                        Igw =   faceSingularityIg(rgj, trit, trit.area, trit.facen̂)*TriGQInfoSglr.weight[gj]
                        Fss    +=   Igw
                        # 对场求积点循环
                    end #for gj
                    Fss /= trit.area
                end
                # Fss项的系数
                Fss *=  lm*abs(arean)

                # 写入数据
                temp = JKη_0div4π*divk²*Fss
                !iszero(δκn)    &&  begin Zmn += δκn*temp;  end
                isbdn           &&  begin Znm += temp;      end
            end
            # 写入数据
            Zts[mi, ni] =   Zmn
            Zst[ni, mi] =   Znm

        end #for mi
    end # for ni 
    return Zts, Zst
end

"""
计算三角形和四面体上相关的 12 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较远的情况
输入
tris    ::  TriangleInfo,  源三角形面体
tetrat  ::  TetrahedraInfo, 场四面体
"""
function EFIEOnNearRWGSWG(tetrat::TetrahedraInfo{IT, FT, CT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    Zst, Zts = EFIEOnNearRWGSWG(tris, tetrat)
    return Zts, Zst
end

# """
# 计算三角形自身相关的 9 个阻抗矩阵元，
# 此函数方法用于计算场源四面体不重合且相隔较远的情况
# 输入
# trit, tetras     :   TetrahedraInfo, 场四面体和源四面体
# """
# function EFIEOnRWGSWG(tri::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
#     Ztt = EFIEOnTris(tri)
#     return Ztt
# end

# """
# 计算四面体自身上相关的 16 个阻抗矩阵元，
# 此函数方法用于计算场源四面体不重合且相隔较远的情况，因此输入有两个四面体信息类型实例
# 输入
# trit, tetras     :   TetrahedraInfo, 场四面体和源四面体
# """
# function EFIEOnRWGSWG(tetrat::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
#     Ztt = EFIEOnTetraSWG(tetrat)
#     return Ztt
# end

# """
# 计算三角形和三角形上相关的 9 个阻抗矩阵元，
# 此函数方法用于计算场源四面体不重合且相隔较远的情况
# 输入
# trit, tetras     :   TetrahedraInfo, 场四面体和源四面体
# """
# function EFIEOnRWGSWG(trit::TriangleInfo{IT, FT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
#     Zts = EFIEOnTris(trit, tris)
#     return Zts, Array(transpose(Zts))
# end

# """
# 计算四面体和四面体上相关的 16 个阻抗矩阵元，
# 此函数方法用于计算场源四面体不重合且相隔较远的情况，因此输入有两个四面体信息类型实例
# 输入
# trit, tetras     :   TetrahedraInfo, 场四面体和源四面体
# """
# function EFIEOnRWGSWG(tetrat::TetrahedraInfo{IT, FT, CT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
#     Zts, Zst = EFIEOnTetrasSWG(tetrat, tetras)
#     return Zts, Zst
# end

# """
# 计算三角形和三角形上相关的 9 个阻抗矩阵元，
# 此函数方法用于计算场源四面体不重合且相隔较近的情况
# 输入
# trit, tetras     :   TetrahedraInfo, 场四面体和源四面体
# """
# function EFIEOnNearRWGSWG(trit::TriangleInfo{IT, FT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
#     Zts = EFIEOnNearTris(trit, tris)
#     return Zts, Array(transpose(Zts))
# end

# """
# 计算四面体和四面体上相关的 16 个阻抗矩阵元，
# 此函数方法用于计算场源四面体不重合且相隔较近的情况，因此输入有两个四面体信息类型实例
# 输入
# trit, tetras     :   TetrahedraInfo, 场四面体和源四面体
# """
# function EFIEOnNearRWGSWG(tetrat::TetrahedraInfo{IT, FT, CT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
#     Zts, Zst = EFIEOnNearTetrasSWG(tetrat, tetras)
#     return Zts, Zst
# end

"""
RWG + SWG 部分的阻抗矩阵
"""
function impedancemat4RWGSWG!(Zmat::Matrix{CT}, trisInfo::AbstractVector{TriangleInfo{IT, FT}}, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 常数
    Rsglr   =   Params.Rsglr
    # 网格元（几何体）数
    trinum      =   length(trisInfo)
    # 线程锁防止对同一数据写入出错
    lockZ   =   SpinLock()
    # 矩阵大小
    nbf     =   size(Zmat, 1)
    # Progress Meter
    pmeter  =   Progress(trinum; desc = "Calculating Z (RWG + SWG) ($nbf × $nbf)...")
    # 外层定义为场基函数循环@threads
    @threads for it in eachindex(trisInfo)
        # 局域的场网格元（几何体）
        @inbounds trit  =   trisInfo[it]
        # 局部判断奇异性距离
        Rsglrlc =   Rsglr
        # 对源网格元（几何体）循环
        @inbounds for js in eachindex(tetrasInfo)
            # 局域的源网格元（几何体）
            tetras    =   tetrasInfo[js]
            # 场源距离
            Rts       =   dist(trit.center, tetras.center)
            # 判断二者远近，调用不同精度的矩阵元处理函数
            if Rts < Rsglrlc
                # 需要进行近奇异性处理的场源网格元（几何体）
                Zts, Zst    =   EFIEOnNearRWGSWG(trit, tetras)
                lock(lockZ)
                for ni in 1:4, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = tetras.inBfsID[ni]
                    # RWG基函数不设半基函数因此跳过
                    ( m == 0 ) && continue
                    # # 往矩阵填充结果
                    Zmat[m, n]  +=  Zts[mi, ni]
                    Zmat[n, m]  +=  Zst[ni, mi]
                end
                unlock(lockZ)
            else
                # 正常高斯求积
                # 计算网格元（几何体）相关的个矩阵元的结果
                Zts, Zst    =   EFIEOnRWGSWG(trit, tetras)
                # 写入数据
                lock(lockZ)
                for ni in 1:4, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = tetras.inBfsID[ni]
                    # RWG基函数不设半基函数因此跳过
                    (m == 0 ) && continue
                    # 往矩阵填充结果
                    Zmat[m, n]  +=  Zts[mi, ni]
                    Zmat[n, m]  +=  Zst[ni, mi]
                end
                unlock(lockZ)

            end # if

        end #for js

        next!(pmeter)

    end #for it
    nothing
end


"""
本函数用于计算金属介质混合体的EFIE阻抗矩阵。
输入信息：
geosInfo    :  为包含几何体信息实例的向量
nbf        :  基函数数目
返回值
Zmat         :  阻抗矩阵

注意，此程序由于采用的对网格元（几何体）循环计算，因此在并行化时，会出现不同线程计算出同一个矩阵元，导致写入冲突，因此要加线程锁保证结果写入正确

"""
function impedancemat4VSIERWGSWG(geosInfo::Vector{VT}, nbf::Integer) where {VT<:AbstractVector}
    
    CT      =   typeof(geosInfo[1][1].ε)
    # 初始化阻抗矩阵
    Zmat    =   zeros(CT, (nbf, nbf))
    # 所有三角形
    trisInfo    =   geosInfo[1]
    # 所有四边形
    tetrasInfo  =   geosInfo[2]
    # PEC + RWG 的部分
    impedancemat4EFIE4PEC!(Zmat, trisInfo, VSBFTypes.sbfType)
    # 介质 + SWG 的部分
    impedancemat4VIE!(Zmat, tetrasInfo, VSBFTypes.vbfType)
    # PEC + 介质 部分
    impedancemat4RWGSWG!(Zmat, trisInfo, tetrasInfo)

    return Zmat
    
end