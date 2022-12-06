
"""
计算三角形上相关9个阻抗矩阵元，
此函数方法用于计算场源三角形不重合且相隔较远的情况，因此输入有两个个三角形信息类型实例
输入
trit， tris     :   TriangleInfo, 场三角形和源三角形
"""
function CFIEOnTris(trit::TriangleInfo{IT, FT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
    CT  =   Complex{FT}

    # 混合场积分方程系数
    α   =   Params.CFIEα
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    Zst     =   zeros(CT, 3, 3)
    # 常数
    JK_0    =   Params.JK_0
    C4divk²     =   Params.C4divk²
    JKηdiv16π   =   Params.JKηdiv16π
    n̂t      =   trit.facen̂
    n̂s      =   tris.facen̂
    # 预分配内存
    gw      =   zero(MMatrix{GQPNTri, GQPNTri, CT})

    # 场源求积点
    rgt     =   getGQPTri(trit)
    rgs     =   getGQPTri(tris)

    ## 求g乘以权重
    # 对源求积点循环
    @inbounds for gj in 1:GQPNTri
        # 源高斯求积点
        rgj  =   view(rgs, :, gj)
        # 对场求积点循环
        for gi in 1:GQPNTri
            # 场高斯求积点
            rgi  =  view(rgt, :, gi)
            # g乘以权重
            gw[gi, gj]  =   greenfunc(rgi, rgj)*TriGQInfo.weight[gi]*TriGQInfo.weight[gj]
        end # for gi
    end #for gj

    # 预分配 ρnj ρmi
    rvec    =   zero(MVec3D{FT})
    ρmi     =   zero(MVec3D{FT})
    ρnj     =   zero(MVec3D{FT})
    # 对源基函数循环求
    @inbounds for ni in 1:3
        # 源基函数序号、带符号边长、自由端
        ln      =   tris.edgel[ni]
        freeVn  =   view(tris.vertices, :, ni)
        # 所在基函数id
        n       =   tris.inBfsID[ni]
        # 判断边是不是基函数（边缘不算）
        n == 0 && continue

        for mi in 1:3
            # 场基函数序号、n带符号边长、自由端
            lm      =   trit.edgel[mi]
            freeVm  =   view(trit.vertices, :, mi)
            # 所在基函数id
            m       =   trit.inBfsID[mi]
            # 判断边是不是基函数（边缘不算）
            m == 0 && continue

            # 用于累加的临时变量
            let ZmnE = zero(CT), ZmnM = zero(CT), ZnmM = zero(CT)
                # 非奇异直接高斯求积
                for gj in 1:GQPNTri
                    # 源高斯求积点
                    rgj     =   view(rgs, :, gj)
                    # ρn
                    ρnj    .=   rgj .- freeVn
                    # 对场求积点循环
                    for gi in 1:GQPNTri
                        # 场高斯求积点
                        rgi     =   view(rgt, :, gi)
                        # ρm
                        ρmi    .=   rgi .- freeVm
                        # rvec
                        rvec   .=   rgi .- rgj
                        # r
                        divr    =   1/norm(rvec)
                        temp    =   (JK_0 + divr) * divr * gw[gi, gj]
                        ρmiρnj  =   ρmi ⋅ ρnj
                        ZmnM   +=   ((ρmi ⋅ rvec) * (n̂t ⋅ ρnj) - (n̂t ⋅ rvec) * ρmiρnj) * temp
                        ZnmM   -=   ((ρnj ⋅ rvec) * (n̂s ⋅ ρmi) - (n̂s ⋅ rvec) * ρmiρnj) * temp
                        # 计算 Zmn\Znm
                        # ZmnM   +=   (ρmi × n̂t) ⋅ (rvec × ρnj) * ((JK_0 + divr) * divr * gw[gi, gj])
                        # ZnmM   -=   (ρnj × n̂s) ⋅ (rvec × ρmi) * ((JK_0 + divr) * divr * gw[gi, gj])                  
                        ZmnE   +=   (ρmi ⋅ ρnj - C4divk²)*gw[gi, gj]
                    end #for gi
                end #for gj
                temp    =   (1-α)*lm*ln*ηdiv16π
                ZmnM   *=   temp
                ZnmM   *=   temp
                ZmnE   *=   α*lm*ln*JKηdiv16π
                Zmn     =   ZmnE + ZmnM
                Znm     =   ZmnE + ZnmM
                # 将结果写入目标数组
                Zts[mi, ni]    =   Zmn
                Zst[ni, mi]    =   Znm
            end # let

        end #for mi
    end # for ni 
    return Zts, Zst
end


"""
计算三角形上相关9个阻抗矩阵元，
此函数方法用于计算场源三角形不重合但相隔较近的情况，因此输入有两个个三角形信息类型实例
输入
trit， tris     :   TriangleInfo, 场三角形和源三角形
"""
function CFIEOnNearTris(trit::TriangleInfo{IT, FT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
    CT  =   Complex{FT}
    # 混合场积分方程系数
    α   =   Params.CFIEα
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    Zst     =   zeros(CT, 3, 3)
    # 常数
    JK_0    =   Params.JK_0
    C4divk²     =   Params.C4divk²
    JKηdiv16π   =   Params.JKηdiv16π
    n̂t      =   trit.facen̂
    n̂s      =   tris.facen̂
    # 预分配内存
    gw      =   zero(MMatrix{GQPNTriSglr, GQPNTriSglr, CT})
    
    # 场源求积点
    rgt     =   getGQPTriSglr(trit)
    rgs     =   getGQPTriSglr(tris)

    ## 求g乘以权重
    # 对源求积点循环
    @inbounds for gj in 1:GQPNTriSglr
        # 源高斯求积点
        rgj  =   view(rgs, :, gj)
        # 对场求积点循环
        for gi in 1:GQPNTriSglr
            # 场高斯求积点
            rgi  =  view(rgt, :, gi)
            # g乘以权重
            gw[gi, gj]  =   greenfunc(rgi, rgj)*TriGQInfoSglr.weight[gi]*TriGQInfoSglr.weight[gj]
        end # for gi
    end #for gj

    # 在场积分点计算得到的 Ig, Ivecg
    Igt     =   zero(MVector{GQPNTriSglr, Complex{FT}})
    Ivecgt  =   zero(MMatrix{3, GQPNTriSglr, Complex{FT}})
    # 循环计算求积点处的奇异项
    for gi in 1:GQPNTriSglr
        Ig, Ivecg   =   faceSingularityIgIvecg(rgt[:, gi], tris, abs(tris.area), tris.facen̂)
        Igt[gi]         =   Ig
        Ivecgt[:, gi]  .=   Ivecg
    end

    # 预分配 ρnj ρmi
    rvec    =   zero(MVec3D{FT})
    ρmi     =   zero(MVec3D{FT})
    ρnj     =   zero(MVec3D{FT})
    # 对源基函数循环求
    @inbounds for ni in 1:3
        # 源基函数序号、带符号边长、自由端
        ln      =   tris.edgel[ni]
        freeVn  =   view(tris.vertices, :, ni)
        # 所在基函数id
        n       =   tris.inBfsID[ni]
        # 判断边是不是基函数（边缘不算）
        n == 0 && continue

        for mi in 1:3
            # 场基函数序号、n带符号边长、自由端
            lm      =   trit.edgel[mi]
            freeVm  =   view(trit.vertices, :, mi)
            # 所在基函数id
            m       =   trit.inBfsID[mi]
            # 判断边是不是基函数（边缘不算）
            m == 0 && continue

            # 用于累加的临时变量
            let ZmnE = zero(CT), ZmnM = zero(CT), ZnmM = zero(CT)
                # 非奇异直接高斯求积
                for gi in 1:GQPNTriSglr
                    # 场高斯求积点
                    rgi     =   view(rgt, :, gi)
                    # ρm
                    ρmi    .=   rgi .- freeVm
                    # ρmin = rgi - freeVn
                    ρmin    =   rgi - freeVn
                    # EFIE
                    ZmnE   +=  ((ρmi ⋅ ρmin - C4divk²)*Igt[gi] - ( ρmi ⋅ view(Ivecgt, :, gi) )) * TriGQInfoSglr.weight[gi]/tris.area
                    # 对源求积点循环
                    for gj in 1:GQPNTriSglr
                        # 源高斯求积点
                        rgj     =   view(rgs, :, gj)
                        # ρn
                        ρnj    .=   rgj .- freeVn
                        # rvec
                        rvec   .=   rgi .- rgj
                        # r
                        divr    =   1/norm(rvec)
                        # 计算 MFIE ZmnM\ZnmM
                        ZmnM   +=   (ρmi × n̂t) ⋅ (rvec × ρnj) * ((JK_0 + divr) * divr * gw[gi, gj])
                        ZnmM   -=   (ρnj × n̂s) ⋅ (rvec × ρmi) * ((JK_0 + divr) * divr * gw[gi, gj])                  
                    end #for gi
                end #for gj
                temp    =   (1-α)*lm*ln*ηdiv16π
                ZmnE   *=   α*lm*ln*JKηdiv16π
                ZmnM   *=   temp
                ZnmM   *=   temp
                Zmn     =   ZmnE + ZmnM
                Znm     =   ZmnE + ZnmM
                # 将结果写入目标数组
                Zts[mi, ni]    =   Zmn
                Zst[ni, mi]    =   Znm
                # Zts[mi, ni]    =   ZmnE
                # Zst[ni, mi]    =   ZmnE
            end # let

        end #for mi
    end # for ni 
    return Zts, Zst
end

"""
计算三角形上相关9个阻抗矩阵元，
此函数方法用于计算场源三角形重合的情况，因此输入只有一个三角形信息类型实例
输入
tri     :   TriangleInfo
"""
function CFIEOnTris(tri::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
    CT  =   Complex{FT}
    # 混合场积分方程系数
    α   =   Params.CFIEα
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)

    # 高斯求积点
    rgt     =   getGQPTriSglr(tri)

    # 常数
    ηdiv8S  =   η_0/(8tri.area)
    # 常数
    C4divk²     =   Params.C4divk²
    JKηdiv16π   =   Params.JKηdiv16π
    # 预分配内存
    gstarw  =   zero(MMatrix{GQPNTriSglr, GQPNTriSglr, Complex{FT}})
    # 预计算奇异项
    F1      =   C4divk²*singularF1(abs(tri.edgel[1]), abs(tri.edgel[2]), abs(tri.edgel[3]))
    
    ## 求gstar乘以权重
    # 对源求积点循环
    @inbounds for gj in 1:GQPNTriSglr
        # 源高斯求积点
        rgj =   view(rgt, :, gj)
        # 对场求积点循环
        for gi in gj:GQPNTriSglr
            # 场高斯求积点
            rgi =   (gi != gj ? view(rgt, :, gi) : rgj)
            # gstar乘以权重
            gstarw[gi, gj]  =   greenfunc_star(rgi, rgj)*TriGQInfoSglr.weight[gi]*TriGQInfoSglr.weight[gj]
        end # for gi
    end #for gj
    
    # 预分配 ρnj ρmi
    rvec    =   zero(MVec3D{FT})
    ρmi     =   zero(MVec3D{FT})
    ρnj     =   zero(MVec3D{FT})
    # 对源基函数循环求
    @inbounds for ni in 1:3
        # 源基函数序号、带符号边长、自由端
        ln      =   tri.edgel[ni]
        freeVn  =   view(tri.vertices, :, ni)
        # 所在基函数id
        n       =   tri.inBfsID[ni]
        # 判断边是不是基函数（边缘不算）
        n == 0 && continue
        
        for mi in ni:3
            # 场基函数序号、n带符号边长、自由端
            lm      =   tri.edgel[mi]
            freeVm  =   view(tri.vertices, :, mi)
            # 所在基函数id
            m       =   tri.inBfsID[mi]
            # 判断边是不是基函数（边缘不算）
            m == 0 && continue

            # 用于累加的临时变量
            ZmnM   =   zero(CT)
            ZmnE   =   zero(CT)
            # 非奇异部分
            for gj in 1:GQPNTriSglr
                # 源高斯求积点
                rgj =   view(rgt, :, gj)
                # ρnj
                ρnj    .=   rgj .- freeVn
                ## EFIE
                # 对场求积点循环
                for gi in gj:GQPNTriSglr
                    # 场高斯求积点
                    rgi =   (gi != gj ? getGQPTriSglr(tri, gi) : rgj)
                    # ρmi
                    ρmi .=   rgi .- freeVm
                    # 计算结果，利用对称性加速计算，因此要对非对角线结果乘2
                    ZmnE   +=  (gi != gj ? 2*(ρmi ⋅ ρnj - C4divk²)*gstarw[gi, gj] : (ρmi ⋅ ρnj - C4divk²)*gstarw[gi, gj])
                end #for gi
                ## MFIE
                # ρnj
                ρmi     .=   rgj .- freeVm
                # 计算
                ZmnM   +=   ρnj ⋅ ρmi * TriGQInfoSglr.weight[gj]
            end #for gj

            # EFIE奇异部分
            if  mi == ni
                # 找出三边
                a, b, c     =   abs.(tri.edgel[MoM_Basics.Vec3IdxCircle[mi:mi+2]])
                #奇异项计算
                ZmnE     +=   singularF21(a, b, c, tri.area^2) - F1
            else
                # 找出三边
                local temp  =   6 - mi - ni
                a, b, c     =   abs.(tri.edgel[MoM_Basics.Vec3IdxCircle[temp:temp+2]])
                # 奇异项计算
                ZmnE     +=   singularF22(a, b, c, tri.area^2) - F1
            end #if

            ZmnM   *=   (1-α)*lm*ln*ηdiv8S
            ZmnE   *=   α*lm*ln*JKηdiv16π
            Zmn     =   ZmnE + ZmnM
            # 将结果写入目标数组
            Zts[mi, ni] =   Zmn
            mi != ni && begin Zts[ni, mi] =   Zmn end
            # Zts[mi, ni]    =   ZmnE 
            # mi != ni && begin Zts[ni, mi]    =   ZmnE end

        end #for mi
    end # for ni 

    return Zts

end

"""本函数用于计算PEC的CFIE阻抗矩阵。
输入信息：
trianglesInfo:  为包含三角形信息实例的向量
nrwg        :  基函数数目
返回值
Zmat         :  阻抗矩阵

注意，此程序由于采用的对三角形循环计算，因此在并行化时，会出现不同线程计算出同一个矩阵元，导致写入冲突，因此要加线程锁保证结果写入正确

"""
function impedancemat4CFIE4PEC(trianglesInfo::Vector{TriangleInfo{IT, FT}}, nrwg::Integer, ::Type{BFT}) where {IT, FT, BFT<:RWG}
    CT = Complex{FT}
    # 初始化阻抗矩阵
    Zmat    =   zeros(CT, (nrwg, nrwg))
    # 索引
    trisIdx =   eachindex(trianglesInfo)
    # 常数
    Rsglr   =   Params.Rsglr
    # 三角形数
    trisnum =   trisIdx.stop
    # 线程锁防止对同一数据写入出错
    lockZ   =   SpinLock()
    # Progress Meter
    pmeter  =   Progress(trisnum, "Calculating Impedance Matrix($nrwg × $nrwg)")

    # 外层定义为场基函数循环
    @threads for triti in trisIdx
        # 局域的场三角形
        @inbounds local trit  =   trianglesInfo[triti]

        @inbounds for trisj in triti:trisnum
            # 局域的源三角形
            local tris  =   trianglesInfo[trisj]
            # 场源距离
            local Rts   =   dist(trit.center, tris.center)

            # 判断二者远近，调用不同精度的矩阵元处理函数
            if Rts == 0.0
                # 计算三角形相关的(3*3)个矩阵元的结果
                Zts  =  CFIEOnTris(trit)
                # 写入数据
                lock(lockZ)
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = trit.inBfsID[ni]
                    # 判断边是不是基函数（边缘不算）
                    (m == 0 || n == 0) && continue
                    # 往矩阵填充结果
                    Zmat[m, n] += Zts[mi, ni]
                end
                unlock(lockZ)

            elseif Rts < Rsglr
                # 需要进行近奇异性处理的场源三角形
                Zts, Zst    =   CFIEOnNearTris(trit, tris)
                
                # 写入数据
                lock(lockZ)
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = tris.inBfsID[ni]

                    # 判断边是不是基函数（边缘不算）
                    (m == 0 || n == 0) && continue
                    
                    Zmat[m, n] += Zts[mi, ni]
                    Zmat[n, m] += Zst[ni, mi]
                end
                unlock(lockZ)

            else
                # 正常高斯求积
                # 计算三角形相关的(3*3)个矩阵元的结果
                Zts, Zst    =   CFIEOnTris(trit, tris)
                
                # 写入数据
                lock(lockZ)
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = tris.inBfsID[ni]

                    # 判断边是不是基函数（边缘不算）
                    (m == 0 || n == 0) && continue
                    
                    Zmat[m, n] += Zts[mi, ni]
                    Zmat[n, m] += Zst[ni, mi]
                end
                unlock(lockZ)

            end # if

        end #for trisj

        next!(pmeter)

    end #for triti

    return Zmat
    
end
