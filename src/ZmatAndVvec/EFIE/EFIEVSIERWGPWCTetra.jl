
# 计算时用到的常数

"""
计算三角形和四面体上相关的 9 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较远的情况
输入
trit    ::  TriangleInfo,  场三角形面体
geos    ::  TetrahedraInfo, 源四面体
"""
function EFIEOnRWGPWC(trit::TriangleInfo{IT, FT}, geos::VT) where {IT<: Integer, FT<:AbstractFloat, VT<:TetrahedraInfo}
    CT  =   Complex{FT} 
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    Zst     =   zeros(CT, 3, 3)
    
    # 场源求积点
    rgt     =   getGQPTri(trit)
    rgs     =   getGQPTetra(geos)

    # 常数项
    JK_0    =   Params.JK_0
    Jη_0divKdVs =   Params.Jη_0divK * geos.volume
    k²      =   Params.k² 

    # 储存并矢结果的临时数组
    dyadG   =   zero(MMatrix{3, 3, CT})
    # 预分配数值以加速
    ρmi     =   zero(MVec3D{FT})
    # 距离向量
    Rtsvec  =   zero(MVec3D{FT})
    # 对源求积点循环
    @inbounds for gi in 1:GQPNTri
        # 场高斯求积点
        rgi  =  view(rgt, :, gi)
        # 并矢赋 0
        dyadG .= 0
        # 对场求积点循环
        for gj in 1:GQPNTetra

            # 源高斯求积点
            rgj  =   view(rgs, :, gj)
            # 距离向量
            Rtsvec .=   rgi .- rgj
            Rts     =   norm(Rtsvec)

            # 1/R
            divR    =   1/Rts
            # 计算 1/R(jk+1/R)
            jkplusR1stdivR1st   =   (JK_0 + divR)*divR

            # R̂R̂并矢
            R̂   =   zero(MVec3D{FT})
            R̂  .=   Rtsvec * divR
            R̂R̂  =   R̂ * R̂'

            # 格林函数项*求积权重
            GreenR  =   exp(-JK_0*Rts)*div4π*divR*TetraGQInfo.weight[gj]
            # 计算矩阵元并叠加
            dyadG .+=     GreenR * ((I - R̂R̂) * k² - (I/3 -  R̂R̂) * (3*jkplusR1stdivR1st))
        end #for gj
        
        for mi in 1:3
            # 所在基函数id
            m       =   trit.inBfsID[mi]
            # 判断边是不是基函数（RWG 边缘不算因此要跳过）
            m  == 0 && continue
            # 场基函数序号、n带符号边长、自由端
            lm      =   trit.edgel[mi]
            freeVm  =   trit.vertices[:, mi]
            # ρm
            ρmi    .=   rgi .- freeVm

            # 累加结果
            temp = (TriGQInfo.weight[gi] * lm /2)
            for ni in 1:3
                Zts[mi, ni] += temp * (ρmi ⋅ dyadG[: , ni])
                Zst[ni, mi] += temp * (ρmi ⋅ dyadG[ni, : ])
            end
        end # mi

    end # for gi
    
    # 常数项修正
    Zts .*= Jη_0divKdVs*geos.κ
    Zst .*= Jη_0divKdVs

    return Zts, Zst
end

"""
计算三角形和四面体上相关的 12 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较远的情况
输入
tris    ::  TriangleInfo,  源三角形面体
geot    ::  TetrahedraInfo, TetrahedraInfo, 场四面体、四面体
"""
function EFIEOnRWGPWC(geot::VT, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat, VT<:VolumeCellType}
    Zst, Zts = EFIEOnRWGPWC(tris, geot)
    return Zts, Zst
end

"""
计算三角形和四面体上相关的 9 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较近的情况
输入
trit    ::  TriangleInfo,  场三角形面体
geos    ::  TetrahedraInfo, TetrahedraInfo, 源四面体、四面体
"""
function EFIEOnNearRWGPWC(trit::TriangleInfo{IT, FT}, geos::VT) where {IT<: Integer, FT<:AbstractFloat, VT<:VolumeCellType}
    CT  =   Complex{FT}
    # 保存结果的小数组
    Zts =   zeros(CT, 3, 3)
    Zst =   zeros(CT, 3, 3)
    # 常数
    Jη_0divK    =   Params.Jη_0divK

    # 场源求积点
    rgt     =   getGQPTriSglr(trit)
    # 预分配数值以加速
    ρmi     =   zero(MVec3D{FT})
    # 储存并矢结果的临时数组
    dyadG   =   zero(MMatrix{3, 3, CT})
    # 对源求积点循环
    @inbounds for gi in 1:GQPNTriSglr
        # 源高斯求积点
        rgi     =   view(rgt, :, gi)
        # 计算 L 算子并矢并乘以权重
        dyadG  .=   volumeSingularityLOpDyad(rgi, geos)
        # 写入结果
        for mi in 1:3
            # 所在基函数id
            m       =   trit.inBfsID[mi]
            # 判断边是不是基函数（RWG 边缘不算因此要跳过）
            m  == 0 && continue
            # 场基函数序号、n带符号边长、自由端
            lm      =   trit.edgel[mi]
            freeVm  =   trit.vertices[:, mi]
            # ρm
            ρmi    .=   rgi .- freeVm

            # 累加结果
            temp = (TriGQInfoSglr.weight[gi] * lm /2)
            for ni in 1:3
                Zts[mi, ni] += temp * (ρmi ⋅ dyadG[: , ni])
                Zst[ni, mi] += temp * (ρmi ⋅ dyadG[ni, : ])
            end
        end # mi
    end

    # 补上常数项
    Zts .*= (Jη_0divK * geos.κ)
    Zst .*= Jη_0divK


    return Zts, Zst
end

"""
计算三角形和四面体上相关的 12 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较远的情况
输入
tris    ::  TriangleInfo,  源三角形面体
geot    ::  TetrahedraInfo, TetrahedraInfo, 场四面体、四面体
"""
function EFIEOnNearRWGPWC(geot::VT, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat, VT<:VolumeCellType}
    Zst, Zts = EFIEOnNearRWGPWC(tris, geot)
    return Zts, Zst
end


"""
RWG + PWC 部分的阻抗矩阵
"""
function impedancemat4RWGPWC!(Zmat::Matrix{CT}, trisInfo::AbstractVector{TriangleInfo{IT, FT}},
    geosInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 常数
    Rsglr   =   Params.Rsglr
    # 网格元（几何体）数
    trinum      =   length(trisInfo)
    # 判断体电流的离散方式，
    discreteJ::Bool = discreteVar == "J"
    # 线程锁防止对同一数据写入出错
    lockZ   =   SpinLock()
    # 矩阵大小
    nbf     =   size(Zmat, 1)
    # Progress Meter
    pmeter  =   Progress(trinum, "Calculating Z (RWG + PWC) ($nbf × $nbf)...")
    # 外层定义为场基函数循环
    @threads for it in eachindex(trisInfo)
        # 局域的场网格元（几何体）@inbounds
        trit  =   trisInfo[it]
        # 局部判断奇异性距离
        Rsglrlc =   Rsglr
        # 对源网格元（几何体）循环@inbounds
        for js in eachindex(geosInfo)
            # 局域的源网格元（几何体）
            geos    =   geosInfo[js]
            κs      =   geos.κ
            # 场源距离
            Rts      =   dist(trit.center, geos.center)
            # 判断二者远近，调用不同精度的矩阵元处理函数
            if Rts < Rsglrlc
                # 需要进行近奇异性处理的场源网格元（几何体）
                Zts, Zst    =   EFIEOnNearRWGPWC(trit, geos)
                lock(lockZ)
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = geos.inBfsID[ni]
                    # RWG基函数不设半基函数因此跳过
                    ( m == 0 ) && continue
                    # # 往矩阵填充结果
                    if discreteJ
                        Zmat[m, n]  +=  Zts[mi, ni]
                    else
                        Zmat[m, n]  +=  Zts[mi, ni] * κs
                    end
                    Zmat[n, m]  +=  Zst[ni, mi]
                end
                unlock(lockZ)
            else
                # 正常高斯求积
                # 计算网格元（几何体）相关的矩阵元的结果
                Zts, Zst    =   EFIEOnRWGPWC(trit, geos)
                # 写入数据
                lock(lockZ)
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = geos.inBfsID[ni]
                    # RWG基函数不设半基函数因此跳过
                    (m == 0 ) && continue
                    # 往矩阵填充结果
                    if discreteJ
                        Zmat[m, n]  +=  Zts[mi, ni]
                    else
                        Zmat[m, n]  +=  Zts[mi, ni] * κs
                    end
                    Zmat[n, m]  +=  Zst[ni, mi]
                end
                unlock(lockZ)

            end # if

        end #for js

        next!(pmeter)

    end #for it
    nothing
end