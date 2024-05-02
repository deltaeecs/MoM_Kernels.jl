
# 计算时用到的常数
using MoM_Basics:getFreeVns, GQ1DID2GQ3DIDVector, getFreeVIDFromGQ3DID

"""
计算六面体上相关的 36 个阻抗矩阵元，
此函数方法用于计算场源六面体不重合且相隔较远的情况，因此输入有两个六面体信息类型实例
输入
hexat, hexas     :   HexahedraInfo, 场六面体和源六面体
"""
function EFIEOnHexasRBF(hexat::HexahedraInfo{IT, FT, CT}, hexas::HexahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts =   zeros(CT, 6, 6)
    Zst =   zeros(CT, 6, 6)
    
    # 场源求积点
    rgt =   getGQPHexa(hexat)
    rgs =   getGQPHexa(hexas)
    # 场源介质对比度
    κt  =   hexat.κ
    κs  =   hexas.κ

    # 场源六面体包含的面
    facest  =   hexat.faces
    facess  =   hexas.faces
    # 常数项
    k²      =   Params.k²
    mJη_0div4πK =   Params.mJη_0div4πK

    ## 求g乘以权重
    # 对源求积点循环
    # 预分配内存
    gw  =   zeros(CT, GQPNHexa, GQPNHexa)
    @inbounds for gj in 1:GQPNHexa
        # 源高斯求积点
        rgj  =   view(rgs, :, gj)
        # 对场求积点循环
        for gi in 1:GQPNHexa
            # 场高斯求积点
            rgi  =  view(rgt, :, gi)
            # g乘以权重
            gw[gi, gj]  =   greenfunc(rgi, rgj)*(HexaGQInfo.weight[gi]*HexaGQInfo.weight[gj])
        end # for gi
    end #for gj

    # F₃ 项与测试、源基函数无关，可以提前计算
    F₃  =   sum(gw)

    # F₄ 项与测试基函数无关，可以提前计算
    F₄s  =   zero(MVector{6, CT})
    # 提前计算该 n 函数所有的自由端
    freeVns =  Vector{MMatrix{3, GQPNQuad, FT, 3*GQPNQuad}}(undef, 6)
    @inbounds for ni in 1:6
        freeVns[ni] =   getFreeVns(hexas, ni)
        # 源基函数带符号面积
        arean   =   hexas.facesArea[ni]
        # 介质对比度变化量
        δκn     =   facess[ni].δκ
        # 是否为半基函数？
        isbdn   =   facess[ni].isbd
        # 第 ni 个面
        face    =   facess[ni]
        if isbdn || ((δκn != 0) && (arean > 0))
            # 保存结果的临时变量
            gtemp   =   zero(CT)
            for gj in 1:GQPNQuad
                # 源高斯求积点
                rgj  =   getGQPQuad(face, gj)
                # 对场求积点循环
                for gi in 1:GQPNHexa
                    # 场高斯求积点
                    rgi  =  rgt[:, gi]
                    # g乘以权重
                    gtemp  +=   greenfunc(rgi, rgj)*(HexaGQInfo.weight[gi]*QuadGQInfo.weight[gj])
                end # for gi
            end #for gj
            # 写入结果
            F₄s[ni] =   gtemp
        end
    end

    # F₅ 项与源基函数无关，可以提前计算
    F₅t  =  zero(MVector{6, CT})
    freeVms =  Vector{MMatrix{3, GQPNQuad, FT, 3*GQPNQuad}}(undef, 6)
    @inbounds for mi in 1:6
        freeVms[mi] =   getFreeVns(hexat, mi)
        # 源基函数带符号面积
        aream   =   hexat.facesArea[mi]
        # 介质对比度变化量
        δκm     =   facest[mi].δκ
        # 是否为半基函数？
        isbdm   =   facest[mi].isbd
        # 第 mi 个面
        face    =   facest[mi]
        if isbdm || ((δκm != 0) && (aream > 0))
            # 保存结果的临时变量
            gtemp   =   zero(CT)
            for gj in 1:GQPNQuad
                # 源高斯求积点
                rgj  =   getGQPQuad(face, gj)
                # 对场求积点循环
                for gi in 1:GQPNHexa
                    # 场高斯求积点
                    rgi  =  view(rgs, :, gi)
                    # g乘以权重
                    gtemp  +=   greenfunc(rgi, rgj)*(HexaGQInfo.weight[gi]*QuadGQInfo.weight[gj])
                end # for gi
            end #for gj
            # 写入结果
            F₅t[mi] =   gtemp
        end
    end

    # 预分配数值以加速
    ρnj     =   zero(MVec3D{FT})
    ρmi     =   zero(MVec3D{FT})
    # 对源基函数循环求
    @inbounds for ni in 1:6
        # 源基函数带符号面积
        arean   =   hexas.facesArea[ni]
        # 介质对比度变化量
        δκn     =   facess[ni].δκ
        # 是否为半基函数？
        isbdn   =   facess[ni].isbd

        # 所在基函数id
        n       =   hexas.inBfsID[ni]
        # F₄项
        F₄      =   F₄s[ni]
        # 提前计算该 n 函数所有的自由端
        freeVnis    =   freeVns[ni]

        for mi in 1:6
            # 场基函数序号、n带符号边长、自由端
            aream   =   hexat.facesArea[mi]
            # 介质对比度变化量
            δκm     =   facest[mi].δκ
            # 是否为半基函数？
            isbdm   =   facest[mi].isbd
            # 所在基函数id
            m       =   hexat.inBfsID[mi]
            # 面积乘积
            aman    =   aream*arean
            # 提前计算该 m 函数所有的自由端
            freeVmis =   freeVms[mi]

            # 初始化阻抗矩阵元
            Zmn     =   zero(CT)
            Znm     =   zero(CT)

            # 系数 C₃
            C₃  =   mJη_0div4πK * aman

            # 计算 F₂项
            F₂  =   zero(CT)
            # 非奇异直接高斯求积
            for gj in 1:GQPNHexa
                # 源高斯求积点
                @views rgj =   rgs[:, gj]
                # 计算源六面体的 “自由端”
                idn3D   =   GQ1DID2GQ3DIDVector[gj]
                idn     =   getFreeVIDFromGQ3DID(idn3D, ni)
                # ρn
                ρnj    .=   rgj .- view(freeVnis, :, idn)
                # 对场求积点循环
                for gi in 1:GQPNHexa
                    # 场高斯求积点
                    @views rgi =   rgt[:, gi]
                    # 计算场六面体的 “自由端”
                    idm3D   =   GQ1DID2GQ3DIDVector[gi]
                    idm     =   getFreeVIDFromGQ3DID(idm3D, mi)
                    # ρn
                    ρmi    .=   rgi .- view(freeVmis, :, idm)
                    # 计算结果
                    F₂  +=  (ρmi ⋅ ρnj) * gw[gi, gj]                    
                end #for gi
            end #for gj

            # CF23, 对称项，Z_mn\Z_nm都用得到
            CF23    =   C₃*(-k²*F₂ + F₃)
            # 累加 Zmn\Znm
            Zmn    +=   κs*CF23
            Znm    +=   κt*CF23

            # F₃ \ F₄ 项已经算出，下面为填充过程
            # F₄ 项，Z[m,n]的 F₄ 项与Z[n,m]的 F₅ 项相等，此处可用来简化计算
            # Zmn-F₄
            (δκn != 0) && (arean > 0) &&    begin Zmn -=    δκn*C₃*F₄; end
            # Znm - F5
            isbdn &&        begin   Znm -=  κt*C₃*F₄;  end

            # F₅ 项，Z[m,n]的 F₅ 项与Z[n,m]的 F₄ 项相等，此处可用来简化计算
            # Zmn-F₅
            isbdm &&        begin   Zmn -=  κs*C₃*F₅t[mi];  end
            # Znm - F₄
            (δκm != 0) && (aream > 0) &&    begin Znm -=    δκm*C₃*F₅t[mi];  end

            # F₆ 项，Z[m,n]的 F₆ 项与Z[n,m]的 F₆ 项相等，此处同样可用来简化计算
            # 判断Z[m,n]的 F₆ 项是否存在的条件
            statem  =   isbdm && (δκn != 0) && (arean>0)
            # 判断Z[n,m]的 F₆ 项是否存在的条件
            staten  =   isbdn && (δκm != 0) && (aream>0)
            if statem || staten
                # 两个面
                facem    =   facest[mi]
                facen    =   facess[ni]
                F₆     =   zero(CT)
                if m != n
                    for gj in 1:GQPNQuad
                        # 源高斯求积点
                        rgj =   getGQPQuad(facen, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNQuad
                            # 场高斯求积点
                            rgi  = getGQPQuad(facem, gi)
                            # g乘以权重
                            F₆  +=   greenfunc(rgi, rgj)*(QuadGQInfo.weight[gi]*QuadGQInfo.weight[gj])
                        end # for gi
                    end #for gj
                else
                    IgSdivS =   zero(CT)
                    # 对场求积点循环
                    for gi in 1:GQPNQuadSglr
                        # 场高斯求积点
                        rgi =   getGQPQuadSglr(facem, gi)
                        # 奇异性Ig
                        Ig  =   faceSingularityIg(rgi, facen, abs(arean), view(hexas.facesn̂, :, ni))
                        # IgSdivS累加
                        IgSdivS +=  QuadGQInfoSglr.weight[gi] * Ig
                    end # for gi
                    # IgSdivS修正并写入结果
                    F₆ +=   IgSdivS/abs(arean)
                end
                # C₃F₆ 项
                C₃F₆    =   C₃*F₆
                # 写入数据
                statem  && begin Zmn += δκn*C₃F₆; end
                staten  && begin Znm += δκm*C₃F₆; end
            end
            # 写入数据
            Zts[mi, ni] =   Zmn
            Zst[ni, mi] =   Znm

        end #for mi
    end # for ni 
    return Zts, Zst
end

"""
计算六面体上相关的 36 个阻抗矩阵元，
此函数方法用于计算场源六面体重合的情况，因此输入有一个六面体信息类型实例
输入
hexat     :   HexahedraInfo, 场六面体和源六面体
"""
function EFIEOnHexaRBF(hexat::HexahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Ztt =   zeros(CT, 6, 6)

    # 源即是场
    hexas  =   hexat
    
    # 场源求积点
    rgt =   getGQPHexaSglr(hexat)
    rgs =   getGQPHexaSSglr(hexat)
    # 场源介质对比度
    κt  =   hexat.κ
    κs  =   κt

    # 常数项
    k²      =   Params.k²
    divJω   =   Params.divJω
    mJη_0div4πK =   Params.mJη_0div4πK

    # divVε
    divVε   =   1/(hexat.volume * hexat.ε)

    # 场源六面体包含的面
    facest  =   hexat.faces
    facess  =   facest

    ## 源体区域的积分 ∫ g(R) dV'
    #  预分配内存
    IgdivVs     =   zero(MVector{GQPNHexaSglr, CT})
    #  对场求积点循环计算上述积分 
    @inbounds for gi in 1:GQPNHexaSglr
        # 场高斯求积点
        rgi =   rgt[:, gi]
        # 对场求积点循环
        IgV =   volumeSingularityIg(rgi, hexas)
        # 写入值
        IgdivVs[gi] =   IgV/hexas.volume
    end #for gi

    # F₃ 项与测试、源基函数无关，可以提前计算
    F₃  =   HexaGQInfoSglr.weight ⋅ IgdivVs

    # F₄ 项与测试基函数无关，可以提前计算
    F₄s =   zero(MVector{6, CT})
    # 自由端也要提前计算，否则会严重拖慢程序
    freeVnsSSglr    =  Vector{MMatrix{3, GQPNQuadSSglr, FT, 3*GQPNQuadSSglr}}(undef, 6)
    freeVns         =  Vector{MMatrix{3, GQPNQuadSglr, FT, 3*GQPNQuadSglr}}(undef, 6)
    freeVms         =  Vector{MMatrix{3, GQPNQuadSglr, FT, 3*GQPNQuadSglr}}(undef, 6)
    @inbounds for ni in 1:6
        freeVnsSSglr[ni]    =   getFreeVnsSSglr(hexas, ni)
        freeVns[ni] =   getFreeVnsSglr(hexas, ni)
        freeVms[ni] =   getFreeVnsSglr(hexat, ni)
        # 该面
        facen   =   facess[ni]
        # 源基函数带符号面积、自由端
        arean   =   hexas.facesArea[ni]
        # 介质对比度变化量
        δκn     =   facen.δκ
        # 是否为半基函数？
        isbdn   =   facen.isbd
        if isbdn || ((δκn != 0) && (arean > 0))
            # 保存结果的临时变量
            IgSdivS =   zero(CT)
            # 对场求积点循环
            for gi in 1:GQPNHexaSglr
                # 场高斯求积点
                rgi =   view(rgt, :, gi)
                # 奇异性Ig
                Ig  =   faceSingularityIg(rgi, facen, abs(arean), view(hexas.facesn̂, :, ni))
                # IgSdivS累加
                IgSdivS +=  HexaGQInfoSglr.weight[gi] * Ig
            end # for gi
            # IgSdivS修正并写入结果
            F₄s[ni] =   IgSdivS/abs(arean)
        end
    end

    # F₅ 项与源基函数无关，可以提前计算，同一体元与 F₄ 相同
    F₅t  =   F₄s
    # 预分配数值以加速
    ρnj     =   zero(MVec3D{FT})
    ρmi     =   zero(MVec3D{FT})
    # 对源基函数循环求积分
    @inbounds for ni in 1:6
        # 源基函数带符号面积、自由端
        arean   =   hexas.facesArea[ni]
        # 介质对比度变化量
        δκn     =   facess[ni].δκ
        # 是否为半基函数？
        isbdn   =   facess[ni].isbd

        # 所在基函数id
        n       =   hexas.inBfsID[ni]
        # F₄项
        F₄      =   F₄s[ni]
        # 提前计算该 n 函数所有的自由端
        freeVnis        =   freeVns[ni]
        freeVnisSSglr   =   freeVnsSSglr[ni]

        for mi in ni:6
            # 场基函数序号、n带符号边长、自由端
            aream   =   hexat.facesArea[mi]
            # 介质对比度变化量
            δκm     =   facest[mi].δκ
            # 是否为半基函数？
            isbdm   =   facest[mi].isbd
            # 所在基函数id
            m       =   hexat.inBfsID[mi]
            # 面积乘积
            aman    =   aream*arean
            # 提前计算该 m 函数所有的自由端
            freeVmis    =   freeVms[mi]

            # 初始化阻抗矩阵元
            Zmn     =   zero(CT)
            Znm     =   zero(CT)

            # C₁
            C₁  =   divJω*divVε * aman
            # 系数 C₃
            C₃  =   mJη_0div4πK * aman

            # 计算 F₁ F₂ 项
            F₁  =   zero(FT)
            F₂  =   zero(CT)
            # 累加
            # 对场点循环
            for gi in 1:GQPNHexaSglr
                # 场高斯求积点
                rgi =   view(rgt, :, gi)
                # 计算场六面体的 “自由端” id
                idm3D   =   GQ1DID2GQ3DIDVectorSglr[gi]
                idm     =   getFreeVIDFromGQ3DIDSglr(idm3D, mi)
                # 计算源六面体的 “自由端” id
                idn     =   getFreeVIDFromGQ3DIDSglr(idm3D, ni)
                # @views freeVn   =   freeVns[:, idn]
                # ρn
                ρnj    .=   rgi .- view(freeVnis, :, idn)
                # ρm
                ρmi    .=   rgi .- view(freeVmis, :, idm)
                # 计算结果
                ρmρn    =   ρmi ⋅ ρnj
                # 计算
                F₁  +=  HexaGQInfoSglr.weight[gi] * ρmρn
                # 对源点循环
                for gj in 1:GQPNHexaSSglr
                    # 计算源六面体的 “自由端” id 
                    idn3D   =   GQ1DID2GQ3DIDVectorSSglr[gj]
                    idn     =   getFreeVIDFromGQ3DIDSSglr(idn3D, mi)
                    # 场高斯求积点
                    rgj     =   view(rgs, :, gj)
                    # ρn
                    ρnj    .=   rgj .- view(freeVnisSSglr, :, idn)
                    # 计算结果
                    ρmρn    =   ρmi ⋅ ρnj
                    # 累加结果
                    F₂  +=  HexaGQInfoSglr.weight[gi] * HexaGQInfoSSglr.weight[gj] * greenfunc(rgi, rgj) * ρmρn
                end #for gj
            end #for gi
            # C₁F₁ 仅在基函数、测试函数定义在同一个六面体时存在
            C₁F₁    =   C₁ * F₁
            # 累加 Zmn\Znm
            Zmn    +=   C₁F₁
            Znm    +=   C₁F₁

            # CF23, 对称项，Z_mn\Z_nm都用得到
            CF23    =   C₃*(-k²*F₂ + F₃)
            # 累加 Zmn\Znm
            Zmn    +=   κs*CF23
            Znm    +=   κt*CF23

            # F₃ \ F₄ 项已经算出，下面为填充过程
            # F₄ 项，Z[m,n]的 F₄ 项与Z[n,m]的 F₅ 项相等，此处可用来简化计算
            # Zmn-F₄
            (δκn != 0) && (arean > 0) &&    begin Zmn   -=  δκn*C₃*F₄s[ni]; end
            # Znm - F5
            isbdn   &&      begin Znm -=    κt*C₃*F₄s[ni];  end

            # F₅ 项，Z[m,n]的 F₅ 项与Z[n,m]的 F₄ 项相等，此处可用来简化计算
            # Zmn-F₅
            isbdm   &&      begin Zmn -=    κs*C₃*F₅t[mi];  end
            # Znm - F₄
            (δκm != 0) && (aream > 0) &&    begin Znm   -=  δκm*C₃*F₅t[mi];  end

            # F₆ 项，Z[m,n]的 F₆ 项与Z[n,m]的 F₆ 项相等，此处同样可用来简化计算
            # 判断Z[m,n]的 F₆ 项是否存在的条件
            statem  =   isbdm && (δκn != 0) && (arean>0)
            # 判断Z[n,m]的 F₆ 项是否存在的条件
            staten  =   isbdn && (δκm != 0) && (aream>0)
            if statem || staten
                # 两个面
                facem    =   facest[mi]
                facen    =   facess[ni]
                F₆     =   zero(CT)
                if m != n
                    for gj in 1:GQPNQuadSglr
                        # 源高斯求积点
                        rgj  =   getGQPQuadSglr(facen, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNQuadSglr
                            # 场高斯求积点
                            rgi  =  getGQPQuadSglr(facem, gi)
                            # g乘以权重
                            F₆  +=   greenfunc(rgi, rgj)*QuadGQInfoSglr.weight[gi]*QuadGQInfoSglr.weight[gj]
                        end # for gi
                    end #for gj
                else
                    IgSdivS =   zero(CT)
                    # 对场求积点循环
                    for gi in 1:GQPNQuadSglr
                        # 场高斯求积点
                        rgi =   getGQPQuadSglr(facem, gi)
                        # 奇异性Ig
                        Ig  =   faceSingularityIg(rgi, facen, abs(arean), view(hexas.facesn̂, :, ni))
                        # IgSdivS累加
                        IgSdivS +=  QuadGQInfoSglr.weight[gi] * Ig
                    end # for gi
                    # IgSdivS修正并写入结果
                    F₆ +=   IgSdivS/abs(arean)
                end
                # C₃F₆ 项
                C₃F₆    =   C₃*F₆
                # 写入数据
                statem  && begin Zmn += δκn*C₃F₆; end
                staten  && begin Znm += δκm*C₃F₆; end

            end
            # 写入数据
            Ztt[mi, ni] =   Zmn
            Ztt[ni, mi] =   Znm

        end #for mi
    end # for ni 
    return Ztt
end

using MoM_Basics:getFreeVnsSglr, getFreeVnsSSglr, GQ1DID2GQ3DIDVectorSglr, GQ1DID2GQ3DIDVectorSSglr
using MoM_Basics:getFreeVIDFromGQ3DIDSglr, getFreeVIDFromGQ3DIDSSglr

"""
计算六面体上相关的 36 个阻抗矩阵元，
此函数方法用于计算场源六面体不重合但相隔较近的情况，输入有两个六面体信息类型实例
输入
hexat, hexas     :   HexahedraInfo, 场六面体和源六面体
"""
function EFIEOnNearHexasRBF(hexat::HexahedraInfo{IT, FT, CT}, hexas::HexahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts =   zeros(CT, 6, 6)
    Zst =   zeros(CT, 6, 6)
    
    # 场源求积点
    rgt =   getGQPHexaSglr(hexat)
    rgs =   getGQPHexaSSglr(hexas)
    # 场源介质对比度
    κt  =   hexat.κ
    κs  =   hexas.κ

    # 场源六面体包含的面
    facest  =   hexat.faces
    facess  =   hexas.faces
    # 常数项
    k²      =   Params.k²
    mJη_0div4πK   =   Params.mJη_0div4πK

    gw  =   zeros(CT, GQPNHexaSglr, GQPNHexaSSglr)
    @inbounds for gj in 1:GQPNHexaSSglr
        # 源高斯求积点
        rgj  =   view(rgs, :, gj)
        # 对场求积点循环
        for gi in 1:GQPNHexaSglr
            # 场高斯求积点
            rgi  =  view(rgt, :, gi)
            # g乘以权重
            gw[gi, gj]  =   greenfunc(rgi, rgj)*(HexaGQInfoSglr.weight[gi]*HexaGQInfoSSglr.weight[gj])
        end # for gi
    end #for gj
    # F₃ 项与测试、源基函数无关，可以提前计算
    F₃  =   sum(gw)

    # F₄ 项与测试基函数无关，可以提前计算
    F₄s =   zero(MVector{6, CT})
    # 自由端也要提前计算，否则会严重拖慢程序
    freeVnsSSglr    =  Vector{MMatrix{3, GQPNQuadSSglr, FT, 3*GQPNQuadSSglr}}(undef, 6)
    freeVms         =  Vector{MMatrix{3, GQPNQuadSglr, FT, 3*GQPNQuadSglr}}(undef, 6)
    @inbounds for ni in 1:6
        freeVnsSSglr[ni]    =   getFreeVnsSSglr(hexas, ni)
        freeVms[ni]         =   getFreeVnsSglr(hexat, ni)
        # 该面
        # 该面
        facen   =   facess[ni]
        # 源基函数带符号面积、自由端
        arean   =   hexas.facesArea[ni]
        # 介质对比度变化量
        δκn     =   facen.δκ
        # 是否为半基函数？
        isbdn   =   facen.isbd
        if isbdn || ((δκn != 0) && (arean > 0))
            # 保存结果的临时变量
            IgSdivS =   zero(CT)
            # 对场求积点循环
            for gi in 1:GQPNHexaSglr
                # 场高斯求积点
                rgi =   rgt[:, gi]
                # 奇异性Ig
                Ig  =   faceSingularityIg(rgi, facen, abs(hexas.facesArea[ni]), view(hexas.facesn̂, :, ni))
                # IgSdivS累加
                IgSdivS +=  HexaGQInfoSglr.weight[gi] * Ig
            end # for gi
            # IgSdivS修正并写入结果
            F₄s[ni] =   IgSdivS/abs(hexas.facesArea[ni])
        end
    end

    # F₅ 项与源基函数无关，可以提前计算
    F₅t  =   zero(MVector{6, CT})
    @inbounds for mi in 1:6
        # 该面
        facem   =   facest[mi]
        # 源基函数带符号面积、自由端
        aream   =   hexat.facesArea[mi]
        # 介质对比度变化量
        δκm     =   facem.δκ
        # 是否为半基函数？
        isbdm   =   facem.isbd
        if isbdm || ((δκm != 0) && (aream > 0))
            # 保存结果的临时变量
            IgSdivS     =   zero(CT)
            # 对场求积点循环
            for gj in 1:GQPNHexaSSglr
                # 场高斯求积点
                rgj  =  view(rgs, :, gj)
                # 奇异性Ig
                Ig  =   faceSingularityIg(rgj, facem, abs(hexat.facesArea[mi]), view(hexat.facesn̂, :, mi))
                # IgSdivS累加
                IgSdivS +=  HexaGQInfoSSglr.weight[gj] * Ig
            end # for gj
            # 写入结果
            F₅t[mi] =   IgSdivS/abs(hexat.facesArea[mi])
        end
    end
    # 预分配数值以加速
    ρnj     =   zero(MVec3D{FT})
    ρmi     =   zero(MVec3D{FT})
    # 对源基函数循环求
    @inbounds for ni in 1:6
        # 源基函数带符号面积、自由端
        arean   =   hexas.facesArea[ni]
        # 介质对比度变化量
        δκn     =   facess[ni].δκ
        # 是否为半基函数？
        isbdn   =   facess[ni].isbd

        # 所在基函数id
        n       =   hexas.inBfsID[ni]
        # F₄项
        F₄      =   F₄s[ni]
        # 提前计算该 n 函数所有的自由端
        freeVnisSSglr   =   freeVnsSSglr[ni]

        for mi in 1:6
            # 场基函数序号、n带符号边长、自由端
            aream   =   hexat.facesArea[mi]
            # 介质对比度变化量
            δκm     =   facest[mi].δκ
            # 是否为半基函数？
            isbdm   =   facest[mi].isbd
            # 所在基函数id
            m       =   hexat.inBfsID[mi]
            # 面积乘积
            aman    =   aream*arean
            # 提前计算该 m 函数所有的自由端
            freeVmis    =   freeVms[mi]

            # 初始化阻抗矩阵元
            Zmn     =   zero(CT)
            Znm     =   zero(CT)

            # 系数 C₃
            C₃  =   mJη_0div4πK * aman

            # 计算 F₂项
            F₂  =   zero(CT)
            # 对场点循环
            for gi in 1:GQPNHexaSglr
                # 场高斯求积点
                rgi =   view(rgt, :, gi)
                # 计算场六面体的 “自由端” id
                idn3D   =   GQ1DID2GQ3DIDVectorSglr[gi]
                idm     =   getFreeVIDFromGQ3DIDSglr(idn3D, mi)
                # ρm
                ρmi    .=   rgi .- view(freeVmis, :, idm)
                # 对源点循环
                for gj in 1:GQPNHexaSSglr
                    # 计算源六面体的 “自由端”
                    # 计算源六面体的 “自由端” id 
                    idn3D   =   GQ1DID2GQ3DIDVectorSSglr[gj]
                    idn     =   getFreeVIDFromGQ3DIDSSglr(idn3D, ni)
                    # 场高斯求积点
                    rgj     =   view(rgs, :, gj)
                    # ρn
                    ρnj    .=   rgj .- view(freeVnisSSglr, :, idn)
                    # 计算结果
                    ρmρn    =   ρmi ⋅ ρnj
                    # 累加结果
                    F₂  +=  HexaGQInfoSglr.weight[gi] * HexaGQInfoSSglr.weight[gj] * greenfunc(rgi, rgj) * ρmρn
                end #for gj
            end #for gi

            # CF23, 对称项，Z_mn\Z_nm都用得到
            CF23    =   C₃*(-k²*F₂ + F₃)
            # 累加 Zmn\Znm
            Zmn    +=   κs*CF23
            Znm    +=   κt*CF23

            # F₃ \ F₄ 项已经算出，下面为填充过程
            # F₄ 项，Z[m,n]的 F₄ 项与Z[n,m]的 F₅ 项相等，此处可用来简化计算
            # Zmn-F₄
            (δκn != 0) && (arean > 0) &&    begin Zmn   -=  δκn*C₃*F₄; end
            # Znm - F5
            isbdn &&        begin Znm -=    κt*C₃*F₄;  end

            # F₅ 项，Z[m,n]的 F₅ 项与Z[n,m]的 F₄ 项相等，此处可用来简化计算
            # Zmn-F₅
            isbdm &&        begin Zmn -=    κs*C₃*F₅t[mi];  end
            # Znm - F₄
            (δκm != 0) && (aream > 0) &&    begin Znm   -=  δκm*C₃*F₅t[mi];  end

            # F₆ 项，Z[m,n]的 F₆ 项与Z[n,m]的 F₆ 项相等，此处同样可用来简化计算
            # 判断Z[m,n]的 F₆ 项是否存在的条件
            statem  =   isbdm && (δκn != 0) && (arean>0)
            # 判断Z[n,m]的 F₆ 项是否存在的条件
            staten  =   isbdn && (δκm != 0) && (aream>0)
            if statem || staten
                # 两个面
                facem    =   facest[mi]
                facen    =   facess[ni]
                F₆     =   zero(CT)
                if m != n
                    for gi in 1:GQPNQuadSglr
                        # 源高斯求积点
                        rgi  =   getGQPQuadSglr(facem, gi)
                        # 对场求积点循环
                        for gj in 1:GQPNQuadSglr
                            # 场高斯求积点
                            rgj  = getGQPQuadSglr(facen, gj)
                            # g乘以权重
                            F₆  +=   greenfunc(rgi, rgj)*QuadGQInfoSglr.weight[gi]*QuadGQInfoSglr.weight[gj]
                        end # for gi
                    end #for gj
                else
                    IgSdivS =   zero(CT)
                    # 对场求积点循环
                    for gi in 1:GQPNQuadSglr
                        # 场高斯求积点
                        rgi =   getGQPQuadSglr(facem, gi)
                        # 奇异性Ig
                        Ig  =   faceSingularityIg(rgi, facen, abs(arean), view(hexas.facesn̂, :, ni))
                        # IgSdivS累加
                        IgSdivS +=  QuadGQInfoSglr.weight[gi] * Ig
                    end # for gi
                    # IgSdivS修正并写入结果
                    F₆ +=   IgSdivS/abs(arean)
                end
                # C₃F₆ 项
                C₃F₆    =   C₃*F₆
                # 写入数据
                statem  && begin Zmn += δκn*C₃F₆; end
                staten  && begin Znm += δκm*C₃F₆; end

            end
            # 写入数据
            Zts[mi, ni] =   Zmn
            Zst[ni, mi] =   Znm

        end #for mi
    end # for ni 
    return Zts, Zst
end

"""
本函数用于计算介质的EFIE阻抗矩阵。
输入信息：
hexasInfo  :  为包含六面体信息实例的向量
nrbf        :  基函数数目
返回值
Zmat         :  阻抗矩阵

注意，此程序由于采用的对六面体循环计算，因此在并行化时，会出现不同线程计算出同一个矩阵元，导致写入冲突，因此要加线程锁保证结果写入正确

"""
function impedancemat4VIE(hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}, nrbf::Integer, bfT::Type{BFT}) where {IT, FT, CT, BFT<:RBF}
    
    # 初始化阻抗矩阵
    Zmat    =   zeros(Complex{FT}, (nrbf, nrbf))
    impedancemat4VIE!(Zmat, hexasInfo, bfT)

    return Zmat
    
end

"""
本函数用于计算介质的EFIE阻抗矩阵。
输入信息：
hexasInfo  :  为包含六面体信息实例的向量
nrbf        :  基函数数目
返回值
Zmat         :  阻抗矩阵

注意，此程序由于采用的对六面体循环计算，因此在并行化时，会出现不同线程计算出同一个矩阵元，导致写入冲突，因此要加线程锁保证结果写入正确

"""
function impedancemat4VIE!(Zmat::Matrix{CT}, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}, ::Type{BFT}) where {IT, FT, CT, BFT<:RBF}
    
    # 常数
    Rsglr       =   Params.Rsglr
    # 六面体数
    hexasnum    =   length(hexasInfo)
    isoffset    =   isa(hexasInfo, OffsetVector)
    geoInterval = begin 
        isoffset ? begin
            st  =   (eachindex(hexasInfo).offset) + 1
            st:(st - 1 + hexasnum)
        end : begin
            1:hexasnum
        end
    end
    # 线程锁防止对同一数据写入出错
    lockZ   =   SpinLock()
    nbf     =   size(Zmat, 1)
    # Progress Meter
    pmeter  =   Progress(hexasnum; desc = "Calculating Z (RBF, EFIE) ($nbf × $nbf)...")
    # 外层定义为场基函数循环
    @threads for it in geoInterval
        # 局域的场六面体
        @inbounds local hexat  =   hexasInfo[it]
        # 局部判断奇异性距离
        Rsglrlc =   Rsglr/sqrt(norm(hexat.ε)/ε_0)
        # 对源六面体循环
        @inbounds for js in it:geoInterval.stop
            # 局域的源六面体
            local hexas    =   hexasInfo[js]
            # 场源距离
            local Rts   =   dist(hexat.center, hexas.center)
            
            if it == js
                # 重合场源六面体
                Ztt     =   EFIEOnHexaRBF(hexat)
                lock(lockZ)
                for ni in 1:6, mi in 1:6
                    # 基函数id
                    m = hexat.inBfsID[mi]
                    n = hexas.inBfsID[ni]
                    # 往矩阵填充结果
                    Zmat[m, n]  +=  Ztt[mi, ni]
                end
                unlock(lockZ)
            # 判断二者远近，调用不同精度的矩阵元处理函数
            elseif Rts < Rsglrlc
                # 需要进行近奇异性处理的场源六面体
                Zts, Zst    =   EFIEOnNearHexasRBF(hexat, hexas)
                lock(lockZ)
                for ni in 1:6, mi in 1:6
                    # 基函数id
                    m = hexat.inBfsID[mi]
                    n = hexas.inBfsID[ni]
                    # 往矩阵填充结果
                    Zmat[m, n]  +=  Zts[mi, ni]
                    Zmat[n, m]  +=  Zst[ni, mi]
                end
                unlock(lockZ)
            else
                # 正常高斯求积
                # 计算六面体相关的(6*6)个矩阵元的结果
                Zts, Zst    =   EFIEOnHexasRBF(hexat, hexas)
                # 写入数据
                lock(lockZ)
                for ni in 1:6, mi in 1:6
                    # 基函数id
                    m = hexat.inBfsID[mi]
                    n = hexas.inBfsID[ni]
                    # 往矩阵填充结果
                    Zmat[m, n]  +=  Zts[mi, ni]
                    Zmat[n, m]  +=  Zst[ni, mi]
                end
                unlock(lockZ)

            end # if

        end #for js

        next!(pmeter)

    end #for it

    return Zmat
    
end