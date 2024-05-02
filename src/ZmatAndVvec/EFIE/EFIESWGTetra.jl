
# 计算时用到的常数

"""
计算四面体上相关的 16 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较远的情况，因此输入有两个四面体信息类型实例
输入
tetrat, tetras     :   TetrahedraInfo, 场四面体和源四面体
"""
function EFIEOnTetrasSWG(tetrat::TetrahedraInfo{IT, FT, CT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts =   zeros(CT, 4, 4)
    Zst =   zeros(CT, 4, 4)
    
    # 场源求积点
    rgt =   getGQPTetra(tetrat)
    rgs =   getGQPTetra(tetras)
    # 场源介质对比度
    κt  =   tetrat.κ
    κs  =   tetras.κ

    # 常数项
    mk²div9     =   Params.mk²div9
    mJη_0div4πK =   Params.mJη_0div4πK

    # 场源四面体包含的面
    facest  =   tetrat.faces
    facess  =   tetras.faces

    ## 求g乘以权重
    # 对源求积点循环
    # 预分配内存
    gw  =   zero(MMatrix{GQPNTetra, GQPNTetra, Complex{FT}})
    @inbounds for gj in 1:GQPNTetra
        # 源高斯求积点
        rgj  =   view(rgs, :, gj)
        # 对场求积点循环
        for gi in 1:GQPNTetra
            # 场高斯求积点
            rgi  =  rgt[:, gi]
            # g乘以权重
            gw[gi, gj]  =   greenfunc(rgi, rgj)*(TetraGQInfo.weight[gi]*TetraGQInfo.weight[gj])
        end # for gi
    end #for gj

    # F₃ 项与测试、源基函数无关，可以提前计算
    F₃  =   sum(gw)

    # F₄ 项与测试基函数无关，可以提前计算
    F₄s  =   zero(MVector{4, CT})
    @inbounds for ni in 1:4
        # 源基函数带符号面积、自由端
        arean   =   tetras.facesArea[ni]
        # 介质对比度变化量
        δκn     =   facess[ni].δκ
        # 是否为半基函数？
        isbdn   =   facess[ni].isbd
        # 第 ni 个面
        faceVert    =   facess[ni]
        if isbdn || ((δκn != 0) && (arean > 0))
            # 保存结果的临时变量
            gtemp   =   zero(CT)
            for gj in 1:GQPNTri
                # 源高斯求积点
                rgj  =   getGQPTri(faceVert, gj)
                # 对场求积点循环
                for gi in 1:GQPNTetra
                    # 场高斯求积点
                    rgi  =  rgt[:, gi]
                    # g乘以权重
                    gtemp  +=   greenfunc(rgi, rgj)*(TetraGQInfo.weight[gi]*TriGQInfo.weight[gj])
                end # for gi
            end #for gj
            # 写入结果
            F₄s[ni] =   gtemp
        end
    end

    # F₅ 项与源基函数无关，可以提前计算
    F₅t  =   zero(MVector{4, CT})
    @inbounds for mi in 1:4
        # 源基函数带符号面积、自由端
        aream   =   tetrat.facesArea[mi]
        # 介质对比度变化量
        δκm     =   facest[mi].δκ
        # 是否为半基函数？
        isbdm   =   facest[mi].isbd
        # 第 mi 个面的三个点
        faceVert    =   facest[mi]
        if isbdm || ((δκm != 0) && (aream > 0))
            # 保存结果的临时变量
            gtemp   =   zero(CT)
            for gj in 1:GQPNTri
                # 源高斯求积点
                rgj  =   getGQPTri(faceVert, gj)
                # 对场求积点循环
                for gi in 1:GQPNTetra
                    # 场高斯求积点
                    rgi  =  view(rgs, :, gi)
                    # g乘以权重
                    gtemp  +=   greenfunc(rgi, rgj)*(TetraGQInfo.weight[gi]*TriGQInfo.weight[gj])
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
    @inbounds for ni in 1:4
        # 源基函数带符号面积、自由端
        arean   =   tetras.facesArea[ni]
        freeVn  =   tetras.vertices[:, ni]
        # 介质对比度变化量
        δκn     =   facess[ni].δκ
        # 是否为半基函数？
        isbdn   =   facess[ni].isbd

        # 所在基函数id
        n       =   tetras.inBfsID[ni]
        # F₄项
        F₄      =   F₄s[ni]

        for mi in 1:4
            # 场基函数序号、n带符号边长、自由端
            aream   =   tetrat.facesArea[mi]
            freeVm  =   tetrat.vertices[:, mi]
            # 介质对比度变化量
            δκm     =   facest[mi].δκ
            # 是否为半基函数？
            isbdm   =   facest[mi].isbd
            # 所在基函数id
            m       =   tetrat.inBfsID[mi]
            # 面积乘积
            aman    =   aream*arean

            # 初始化阻抗矩阵元
            Zmn     =   zero(CT)
            Znm     =   zero(CT)

            # 系数 C₃
            C₃  =   mJη_0div4πK * aman

            # 计算 F₂项
            F₂  =   zero(CT)
            # 非奇异直接高斯求积
            for gj in 1:GQPNTetra
                # 源高斯求积点
                @views rgj =   rgs[:, gj]
                # ρn
                ρnj    .=   rgj .- freeVn
                # if (mi==2) & (ni==1)
                #     println("ρnj ", ρnj, "\t rgj", rgj, "\t freevn", freeVn)
                # end
                # 对场求积点循环
                for gi in 1:GQPNTetra
                    # 场高斯求积点
                    # rgi =  rgt[:, gi]
                    # # ρm
                    # ρmi =  rgi - freeVm
                    @views rgi =   rgt[:, gi]
                    # ρn
                    ρmi    .=   rgi .- freeVm
                    # 计算结果
                    F₂  +=  (ρmi ⋅ ρnj) * gw[gi, gj]                    
                end #for gi
            end #for gj

            # CF23, 对称项，Z_mn\Z_nm都用得到
            CF23    =   C₃*(mk²div9*F₂ + F₃)
            # 累加 Zmn\Znm
            Zmn    +=   κs*CF23
            Znm    +=   κt*CF23

            # F₃ \ F₄ 项已经算出，下面为填充过程
            # F₄ 项，Z[m,n]的 F₄ 项与Z[n,m]的 F₅ 项相等，此处可用来简化计算
            # Zmn-F₄
            (δκn != 0) && (arean > 0) &&    begin Zmn -=    δκn*C₃*F₄s[ni]; end
            # Znm - F5
            isbdn &&        begin   Znm -=  κt*C₃*F₄s[ni];  end

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
                # 两个面的角点
                faceVertt    =   facest[mi]
                faceVerts    =   facess[ni]
                F₆     =   zero(CT)
                if m != n
                    for gj in 1:GQPNTri
                        # 源高斯求积点
                        rgj  =   getGQPTri(faceVertt, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNTri
                            # 场高斯求积点
                            rgi  = getGQPTri(faceVerts, gi)
                            # g乘以权重
                            F₆  +=   greenfunc(rgi, rgj)*(TriGQInfo.weight[gi]*TriGQInfo.weight[gj])
                        end # for gi
                    end #for gj
                else
                    for gj in 1:GQPNTriSglr
                        # 源高斯求积点
                        rgj  =      getGQPTriSglr(faceVertt, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNTriSglr
                            # 场高斯求积点
                            rgi  = getGQPTriSglr(faceVerts, gi)
                            # g乘以权重
                            F₆  +=   greenfunc_star(rgi, rgj)*(TriGQInfoSglr.weight[gi]*TriGQInfoSglr.weight[gj])
                        end # for gi
                    end #for gj
                    F₆  +=     singularF1(facess[ni].edgel...)
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
计算四面体上相关的 16 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较远的情况，因此输入有两个四面体信息类型实例
输入
tetrat, tetras     :   TetrahedraInfo, 场四面体和源四面体
"""
function EFIEOnTetraSWG(tetrat::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Ztt =   zeros(CT, 4, 4)

    # 源即是场
    tetras  =   tetrat
    
    # 场源求积点
    rgt =   getGQPTetraSglr(tetrat)
    rgs =   rgt
    # 场源介质对比度
    κt  =   tetrat.κ
    κs  =   κt
    
    # 常数项
    mk²div9     =   Params.mk²div9
    mJη_0div4πK =   Params.mJη_0div4πK
    div9Jω      =   Params.div9Jω
    # divVε
    divVε   =   1/(tetrat.volume * tetrat.ε)

    # 场源四面体包含的面
    facest  =   tetrat.faces
    facess  =   facest

    ## 源体区域的积分 ∫ g(R) dV' ∫ Rvec g(R) dV'
    #  预分配内存
    IgdivVs     =   zero(MVector{GQPNTetraSglr, CT})
    IvecgdivVs  =   zero(MMatrix{3, GQPNTetraSglr, CT})
    #  对源求积点循环计算上述积分 
    @inbounds for gi in 1:GQPNTetraSglr
        # 源高斯求积点
        rgi  =   rgt[:, gi]
        # 对场求积点循环
        IgV, IvecgV =   volumeSingularityIgIvecg(rgi, tetras)
        # 写入值
        IgdivVs[gi]         =   IgV/tetras.volume
        IvecgdivVs[:, gi]  .=   IvecgV/tetras.volume
    end #for gj

    # F₃ 项与测试、源基函数无关，可以提前计算
    F₃  =   TetraGQInfoSglr.weight ⋅ IgdivVs

    # F₄ 项与测试基函数无关，可以提前计算
    F₄s =   zero(MVector{4, CT})
    @inbounds for ni in 1:4
        # 该面
        facen   =   facess[ni]
        # 源基函数带符号面积、自由端
        arean   =   tetras.facesArea[ni]
        # 介质对比度变化量
        δκn     =   facen.δκ
        # 是否为半基函数？
        isbdn   =   facen.isbd
        if isbdn || ((δκn != 0) && (arean > 0))
            # 保存结果的临时变量
            IgSdivS =   zero(CT)
            # 对场求积点循环
            for gi in 1:GQPNTetraSglr
                # 场高斯求积点
                rgi =   view(rgt, :, gi)
                # 奇异性Ig
                Ig  =   faceSingularityIg(rgi, facen, abs(arean), view(tetras.facesn̂, :, ni))
                # IgSdivS累加
                IgSdivS +=  TetraGQInfoSglr.weight[gi] * Ig
            end # for gi
            # IgSdivS修正并写入结果
            F₄s[ni] =   IgSdivS/abs(arean)
        end
    end

    # F₅ 项与源基函数无关，可以提前计算，同一体元与 F₄ 相同
    F₅t  =   F₄s
        
    # 对源基函数循环求
    @inbounds for ni in 1:4
        # 源基函数带符号面积、自由端
        arean   =   tetras.facesArea[ni]
        freeVn  =   tetras.vertices[:, ni]
        # 介质对比度变化量
        δκn     =   facess[ni].δκ
        # 是否为半基函数？
        isbdn   =   facess[ni].isbd

        # 所在基函数id
        n       =   tetras.inBfsID[ni]
        # F₄项
        F₄      =   F₄s[ni]

        for mi in ni:4
            # 场基函数序号、n带符号边长、自由端
            aream   =   tetrat.facesArea[mi]
            freeVm  =   tetrat.vertices[:, mi]
            # 介质对比度变化量
            δκm     =   facest[mi].δκ
            # 是否为半基函数？
            isbdm   =   facest[mi].isbd
            # 所在基函数id
            m       =   tetrat.inBfsID[mi]
            # 面积乘积
            aman    =   aream*arean

            # 初始化阻抗矩阵元
            Zmn     =   zero(CT)
            Znm     =   zero(CT)

            # C₁
            C₁  =   div9Jω*divVε*aman
            # 系数 C₃
            C₃  =   mJη_0div4πK * aman

            # 计算 F₁ F₂ 项
            F₁  =   zero(FT)
            F₂  =   zero(CT)
            # 累加
            for gi in 1:GQPNTetraSglr
                # 场高斯求积点
                rgi =   rgt[:, gi]
                # ρm
                ρmi =   rgi - freeVm
                # ρmin = rgi - freeVn
                ρmin=   rgi - freeVn
                # 计算结果
                ρmρn=   ρmi ⋅ ρmin
                F₁  +=  TetraGQInfoSglr.weight[gi] * ρmρn
                F₂  +=  TetraGQInfoSglr.weight[gi] * (- ρmi ⋅ view(IvecgdivVs, :, gi) + ρmρn * IgdivVs[gi])
            end #for gi
            # C₁F₁ 仅在基函数、测试函数定义在同一个四面体时存在
            C₁F₁    =   C₁ * F₁
            # 累加 Zmn\Znm
            Zmn    +=   C₁F₁
            Znm    +=   C₁F₁

            # CF23, 对称项，Z_mn\Z_nm都用得到
            CF23    =   C₃*(mk²div9*F₂ + F₃)
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
                # 两个面的角点
                faceVertt    =   facest[mi]
                faceVerts    =   facess[ni]
                F₆     =   zero(CT)
                if m != n
                    for gj in 1:GQPNTriSglr
                        # 源高斯求积点
                        rgj  =   getGQPTriSglr(faceVertt, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNTriSglr
                            # 场高斯求积点
                            rgi  = getGQPTriSglr(faceVerts, gi)
                            # g乘以权重
                            F₆  +=   greenfunc(rgi, rgj)*TriGQInfoSglr.weight[gi]*TriGQInfoSglr.weight[gj]
                        end # for gi
                    end #for gj
                else
                    for gj in 1:GQPNTriSglr
                        # 源高斯求积点
                        rgj  =      getGQPTriSglr(faceVertt, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNTriSglr
                            # 场高斯求积点
                            rgi  =  getGQPTriSglr(faceVerts, gi)
                            # g乘以权重
                            F₆  +=  greenfunc_star(rgi, rgj)*TriGQInfoSglr.weight[gi]*TriGQInfoSglr.weight[gj]
                        end # for gi
                    end #for gj
                    F₆  +=     singularF1(facess[ni].edgel...)  
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

"""
计算四面体上相关的 16 个阻抗矩阵元，
此函数方法用于计算场源四面体不重合但相隔较近的情况，输入有两个四面体信息类型实例
输入
tetrat, tetras     :   TetrahedraInfo, 场四面体和源四面体
"""
function EFIEOnNearTetrasSWG(tetrat::TetrahedraInfo{IT, FT, CT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts =   zeros(CT, 4, 4)
    Zst =   zeros(CT, 4, 4)
    
    # 场源求积点
    rgt =   getGQPTetraSglr(tetrat)
    rgs =   getGQPTetraSglr(tetras)
    # 场源介质对比度
    κt  =   tetrat.κ
    κs  =   tetras.κ
    # 常数项
    mk²div9     =   Params.mk²div9
    mJη_0div4πK =   Params.mJη_0div4πK

    # 场源四面体包含的面
    facest  =   tetrat.faces
    facess  =   tetras.faces

    ## 源体区域的积分 ∫ g(R) dV' ∫ Rvec g(R) dV'
    #  预分配内存
    IgdivVs     =   zero(MVector{GQPNTetraSglr, CT})
    IvecgdivVs  =   zero(MMatrix{3, GQPNTetraSglr, CT})
    #  对源求积点循环计算上述积分 
    @inbounds for gi in 1:GQPNTetraSglr
        # 源高斯求积点
        rgi  =   rgt[:, gi]
        # 对场求积点循环
        IgV, IvecgV =   volumeSingularityIgIvecg(rgi, tetras)
        # 写入值
        IgdivVs[gi]         =   IgV/tetras.volume
        IvecgdivVs[:, gi]  .=   IvecgV/tetras.volume
    end #for gj

    # F₃ 项与测试、源基函数无关，可以提前计算
    F₃  =   TetraGQInfoSglr.weight ⋅ IgdivVs

    # F₄ 项与测试基函数无关，可以提前计算
    F₄s =   zero(MVector{4, CT})
    @inbounds for ni in 1:4
        # 该面
        facen   =   facess[ni]
        # 源基函数带符号面积、自由端
        arean   =   tetras.facesArea[ni]
        # 介质对比度变化量
        δκn     =   facen.δκ
        # 是否为半基函数？
        isbdn   =   facen.isbd
        if isbdn || ((δκn != 0) && (arean > 0))
            # 保存结果的临时变量
            IgSdivS =   zero(CT)
            # 对场求积点循环
            for gi in 1:GQPNTetraSglr
                # 场高斯求积点
                rgi =   rgt[:, gi]
                # 奇异性Ig
                Ig  =   faceSingularityIg(rgi, facen, abs(tetras.facesArea[ni]), view(tetras.facesn̂, :, ni))
                # IgSdivS累加
                IgSdivS +=  TetraGQInfoSglr.weight[gi] * Ig
            end # for gi
            # IgSdivS修正并写入结果
            F₄s[ni] =   IgSdivS/abs(tetras.facesArea[ni])
        end
    end

    # F₅ 项与源基函数无关，可以提前计算
    F₅t  =   zero(MVector{4, CT})
    @inbounds for mi in 1:4
        # 该面
        facem   =   facest[mi]
        # 源基函数带符号面积、自由端
        aream   =   tetrat.facesArea[mi]
        # 介质对比度变化量
        δκm     =   facem.δκ
        # 是否为半基函数？
        isbdm   =   facem.isbd
        if isbdm || ((δκm != 0) && (aream > 0))
            # 保存结果的临时变量
            IgSdivS     =   zero(CT)
            # 对场求积点循环
            for gj in 1:GQPNTetraSglr
                # 场高斯求积点
                rgj  =  view(rgs, :, gj)
                # 奇异性Ig
                Ig  =   faceSingularityIg(rgj, facem, abs(tetrat.facesArea[mi]), view(tetrat.facesn̂, :, mi))
                # IgSdivS累加
                IgSdivS +=  TetraGQInfoSglr.weight[gj] * Ig
            end # for gj
            # 写入结果
            F₅t[mi] =   IgSdivS/abs(tetrat.facesArea[mi])
        end
    end

    # 对源基函数循环求
    @inbounds for ni in 1:4
        # 源基函数带符号面积、自由端
        arean   =   tetras.facesArea[ni]
        freeVn  =   tetras.vertices[:, ni]
        # 介质对比度变化量
        δκn     =   facess[ni].δκ
        # 是否为半基函数？
        isbdn   =   facess[ni].isbd

        # 所在基函数id
        n       =   tetras.inBfsID[ni]
        # F₄项
        F₄      =   F₄s[ni]

        for mi in 1:4
            # 场基函数序号、n带符号边长、自由端
            aream   =   tetrat.facesArea[mi]
            freeVm  =   tetrat.vertices[:, mi]
            # 介质对比度变化量
            δκm     =   facest[mi].δκ
            # 是否为半基函数？
            isbdm   =   facest[mi].isbd
            # 所在基函数id
            m       =   tetrat.inBfsID[mi]
            # 面积乘积
            aman    =   aream*arean

            # 初始化阻抗矩阵元
            Zmn     =   zero(CT)
            Znm     =   zero(CT)

            # 系数 C₃
            C₃  =   mJη_0div4πK * aman

            # 计算 F₂项
            F₂  =   zero(CT)
            # 奇异项累加
            for gi in 1:GQPNTetraSglr
                # 场高斯求积点
                rgi =   rgt[:, gi]
                # ρm
                ρmi =   rgi - freeVm
                # ρmin = rgi - freeVn
                ρmin=   rgi - freeVn
                # 计算结果
                F₂  +=  TetraGQInfoSglr.weight[gi] * (- ρmi ⋅ view(IvecgdivVs, :, gi) + ρmi ⋅ ρmin * IgdivVs[gi])
            end #for gi

            # CF23, 对称项，Z_mn\Z_nm都用得到
            CF23    =   C₃*(mk²div9*F₂ + F₃)
            # 累加 Zmn\Znm
            Zmn    +=   κs*CF23
            Znm    +=   κt*CF23

            # F₃ \ F₄ 项已经算出，下面为填充过程
            # F₄ 项，Z[m,n]的 F₄ 项与Z[n,m]的 F₅ 项相等，此处可用来简化计算
            # Zmn-F₄
            (δκn != 0) && (arean > 0) &&    begin Zmn   -=  δκn*C₃*F₄s[ni]; end
            # Znm - F5
            isbdn &&        begin Znm -=    κt*C₃*F₄s[ni];  end

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
                # 两个面的角点
                faceVertt    =   facest[mi]
                faceVerts    =   facess[ni]
                F₆     =   zero(CT)
                if m != n
                    for gj in 1:GQPNTriSglr
                        # 源高斯求积点
                        rgj  =   getGQPTriSglr(faceVertt, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNTriSglr
                            # 场高斯求积点
                            rgi  = getGQPTriSglr(faceVerts, gi)
                            # g乘以权重
                            F₆  +=   greenfunc(rgi, rgj)*TriGQInfoSglr.weight[gi]*TriGQInfoSglr.weight[gj]
                        end # for gi
                    end #for gj
                else
                    for gj in 1:GQPNTriSglr
                        # 源高斯求积点
                        rgj  =      getGQPTriSglr(faceVertt, gj)
                        # 对场求积点循环
                        for gi in 1:GQPNTriSglr
                            # 场高斯求积点
                            rgi  = getGQPTriSglr(faceVerts, gi)
                            # g乘以权重
                            F₆  +=   greenfunc_star(rgi, rgj)*TriGQInfoSglr.weight[gi]*TriGQInfoSglr.weight[gj]
                        end # for gi
                    end #for gj
                    F₆  +=     singularF1(facess[ni].edgel...)  
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
tetrasInfo  :  为包含四面体信息实例的向量
nswg        :  基函数数目
返回值
Zmat         :  阻抗矩阵

注意，此程序由于采用的对四面体循环计算，因此在并行化时，会出现不同线程计算出同一个矩阵元，导致写入冲突，因此要加线程锁保证结果写入正确

"""
function impedancemat4VIE(tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, nswg::Integer, bfT::Type{BFT}) where {IT, FT, CT, BFT<:SWG}
    
    # 初始化阻抗矩阵
    Zmat    =   zeros(Complex{FT}, (nswg, nswg))
    # 计算
    impedancemat4VIE!(Zmat, tetrasInfo, bfT)

    return Zmat
    
end

"""
本函数用于计算介质的 EFIE 阻抗矩阵。
输入信息：
Zmat       :   已初始化的阻抗矩阵
tetrasInfo  :  为包含四面体信息实例的向量
nswg        :  基函数数目
返回值
Zmat         :  阻抗矩阵

注意，此程序由于采用的对四面体循环计算，因此在并行化时，会出现不同线程计算出同一个矩阵元，导致写入冲突，因此要加线程锁保证结果写入正确

"""
function impedancemat4VIE!(Zmat::Matrix{CT}, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, ::Type{BFT}) where {IT, FT, CT, BFT<:SWG}
    
    # 索引
    tetrasIdx   =   eachindex(tetrasInfo)
    # 常数
    Rsglr       =   Params.Rsglr
    # 四面体数
    tetrasnum   =   length(tetrasInfo)
    # 处理混合网格
    isoffset    =   isa(tetrasInfo, OffsetVector)
    geoInterval = begin 
        isoffset ? begin
                st  =   (eachindex(tetrasInfo).offset) + 1
                st:(st - 1 + tetrasnum)
            end : begin
            1:tetrasnum
        end
    end

    # 线程锁防止对同一数据写入出错
    lockZ   =   SpinLock()
    # 矩阵大小
    nbf     =   size(Zmat, 1)
    # Progress Meter
    pmeter  =   Progress(tetrasnum; desc = "Calculating Z (SWG, EFIE) ($nbf × $nbf)")
    # 外层定义为场基函数循环
    @threads for it in geoInterval
        # 局域的场四面体
        @inbounds local tetrat  =   tetrasInfo[it]
        # 局部判断奇异性距离
        Rsglrlc =   Rsglr/sqrt(norm(tetrat.ε)/ε_0)
        # 对源四面体循环
        @inbounds for js in it:geoInterval.stop
            # 局域的源四面体
            local tetras    =   tetrasInfo[js]
            # 场源距离
            local Rts   =   dist(tetrat.center, tetras.center)
            # isapprox(Rts, Rsglrlc, rtol = 1e-2) && @show it, js
            
            if it == js
                # 重合场源四面体
                Ztt     =   EFIEOnTetraSWG(tetrat)
                lock(lockZ)
                for ni in 1:4, mi in 1:4
                    # 基函数id
                    m = tetrat.inBfsID[mi]
                    n = tetras.inBfsID[ni]
                    # 往矩阵填充结果
                    Zmat[m, n]  +=  Ztt[mi, ni]
                end
                unlock(lockZ)
            # 判断二者远近，调用不同精度的矩阵元处理函数
            elseif Rts < Rsglrlc
                # 需要进行近奇异性处理的场源四面体
                Zts, Zst    =   EFIEOnNearTetrasSWG(tetrat, tetras)
                lock(lockZ)
                for ni in 1:4, mi in 1:4
                    # 基函数id
                    m = tetrat.inBfsID[mi]
                    n = tetras.inBfsID[ni]
                    # 往矩阵填充结果
                    Zmat[m, n]  +=  Zts[mi, ni]
                    Zmat[n, m]  +=  Zst[ni, mi]
                end
                unlock(lockZ)
            else
                # 正常高斯求积
                # 计算四面体相关的(4*4)个矩阵元的结果
                Zts, Zst    =   EFIEOnTetrasSWG(tetrat, tetras)
                # 写入数据
                lock(lockZ)
                for ni in 1:4, mi in 1:4
                    # 基函数id
                    m = tetrat.inBfsID[mi]
                    n = tetras.inBfsID[ni]
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