
"""
计算三角形上相关9个阻抗矩阵元，
此函数方法用于计算场源三角形不重合且相隔较远的情况，因此输入有两个个三角形信息类型实例
输入
trit， tris     :   TriangleInfo, 场三角形和源三角形
"""
function EFIEOnTris(trit::TriangleInfo{IT, FT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
    # 保存结果的 (3*3) 小数组
    Ztri    =   zeros(Complex{FT}, 3, 3)
    # 常数
    C4divk²     =   Params.C4divk²
    JKηdiv16π   =   Params.JKηdiv16π
    # 预分配内存
    gw      =   zero(MMatrix{GQPNTri, GQPNTri, Complex{FT}})
    
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
        
    # 对源基函数循环求
    @inbounds for ni in 1:3
        # 源基函数序号、带符号边长、自由端
        ln      =   tris.edgel[ni]
        freeVn  =   tris.vertices[:, ni]
        # 所在基函数id
        n       =   tris.inBfsID[ni]
        # 判断边是不是基函数（边缘不算）
        n == 0 && continue

        for mi in 1:3
            # 场基函数序号、n带符号边长、自由端
            lm      =   trit.edgel[mi]
            freeVm  =   trit.vertices[:, mi]
            # 所在基函数id
            m       =   trit.inBfsID[mi]
            # 判断边是不是基函数（边缘不算）
            m == 0 && continue

            # 用于累加的临时变量
            let Ztemp::Complex{FT} = 0im
                # 非奇异直接高斯求积
                for gj in 1:GQPNTri
                    # 源高斯求积点
                    rgj =   rgs[:, gj]
                    # ρn
                    ρnj  =   rgj - freeVn
                    # if (mi==2) & (ni==1)
                    #     println("ρnj ", ρnj, "\t rgj", rgj, "\t freevn", freeVn)
                    # end
                    # 对场求积点循环
                    for gi in 1:GQPNTri
                        # 场高斯求积点
                        rgi =  rgt[:, gi]
                        # ρm
                        ρmi =  rgi - freeVm
                        # 计算结果
                        Ztemp  +=  (ρmi ⋅ ρnj - C4divk²)*gw[gi, gj]                    
                    end #for gi
                end #for gj
                Ztemp   *=  lm*ln*JKηdiv16π
                # 将结果写入目标数组
                Ztri[mi, ni]    =   Ztemp
            end # let

        end #for mi
    end # for ni 
    return Ztri
end

"""
计算三角形上相关9个阻抗矩阵元，
此函数方法用于计算场源三角形不重合但相隔较近的情况，因此输入有两个个三角形信息类型实例
输入
trit， tris     :   TriangleInfo, 场三角形和源三角形
"""
function EFIEOnNearTris(trit::TriangleInfo{IT, FT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
    # 保存结果的 (3*3) 小数组
    Ztri    =   zeros(Complex{FT}, 3, 3)
    # 常数
    C4divk²     =   Params.C4divk²
    JKηdiv16π   =   Params.JKηdiv16π
    # 场三角形的高斯求积点
    rgt     =   getGQPTriSglr(trit)
    # 在场积分点计算得到的 Ig, Ivecg
    Igt     =   zero(MVector{GQPNTriSglr, Complex{FT}})
    Ivecgt  =   zero(MMatrix{3, GQPNTriSglr, Complex{FT}})
    # 循环计算求积点处的奇异项
    for gi in 1:GQPNTriSglr
        Ig, Ivecg   =   faceSingularityIgIvecg(rgt[:, gi], tris, abs(tris.area), tris.facen̂)
        Igt[gi]         =   Ig
        Ivecgt[:, gi]  .=   Ivecg
    end
    
    # 对源基函数循环求阻抗矩阵元
    @inbounds for ni in 1:3
        # 源基函数序号、带符号边长、自由端
        ln      =   tris.edgel[ni]
        freeVn  =   tris.vertices[:, ni]
        # 所在基函数id
        n       =   tris.inBfsID[ni]
        # 判断边是不是基函数（边缘半基函数不算）
        n == 0 && continue

        for mi in 1:3
            # 场基函数序号、n带符号边长、自由端
            lm      =   trit.edgel[mi]
            freeVm  =   trit.vertices[:, mi]
            # 所在基函数id
            m       =   trit.inBfsID[mi]
            # 判断边是不是基函数（边缘半基函数不算）
            m == 0 && continue

            # 用于累加的临时变量
            let Ztemp = zero(Complex{FT})
                # 奇异部分
                # 对场求积点循环
                for gi in 1:GQPNTriSglr
                    # 场高斯求积点
                    rgi =  rgt[:, gi]
                    # ρm
                    ρmi =   rgi - freeVm
                    # ρmin = rgi - freeVn
                    ρmin=   rgi - freeVn
                    Ztemp   +=  ((ρmi ⋅ ρmin - C4divk²)*Igt[gi] - ( ρmi ⋅ view(Ivecgt, :, gi) )) * TriGQInfoSglr.weight[gi]/tris.area
                end
                Ztemp   *=  lm*ln*JKηdiv16π
                # 将结果写入目标数组
                Ztri[mi, ni]    =   Ztemp
            end # let

        end #for mi
    end # for ni 
    return Ztri
end

"""
计算三角形上相关9个阻抗矩阵元，
此函数方法用于计算场源三角形重合的情况，因此输入只有一个三角形信息类型实例
输入
tri     :   TriangleInfo
"""
function EFIEOnTris(tri::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
    CT  =   Complex{FT}
    # 面积平方
    local areasquare    =  tri.area^2
    # 保存结果的 (3*3) 小数组
    Ztri    =   zeros(Complex{FT}, 3, 3)
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
        rgj =   getGQPTriSglr(tri, gj)
        # 对场求积点循环
        for gi in gj:GQPNTriSglr
            # 场高斯求积点
            rgi =   (gi != gj ? getGQPTriSglr(tri, gi) : rgj)
            # gstar乘以权重
            gstarw[gi, gj]  =   greenfunc_star(rgi, rgj)*TriGQInfoSglr.weight[gi]*TriGQInfoSglr.weight[gj]
        end # for gi
    end #for gj
        
    # 对源基函数循环求
    @inbounds for ni in 1:3
        # 源基函数序号、带符号边长、自由端
        ln      =   tri.edgel[ni]
        freeVn  =   tri.vertices[:, ni]
        # 所在基函数id
        n       =   tri.inBfsID[ni]
        # 判断边是不是基函数（边缘不算）
        n == 0 && continue
        for mi in ni:3
            # 场基函数序号、n带符号边长、自由端
            lm      =   tri.edgel[mi]
            freeVm  =   tri.vertices[:, mi]
            # 所在基函数id
            m       =   tri.inBfsID[mi]
            # 判断边是不是基函数（边缘不算）
            m == 0 && continue

            # 用于累加的临时变量
            Ztemp = zero(CT)
            # 非奇异部分
            for gj in 1:GQPNTriSglr
                # 源高斯求积点
                rgj =   getGQPTriSglr(tri, gj)
                # ρs
                ρs  =   rgj - freeVn
                # 对场求积点循环
                for gi in gj:GQPNTriSglr
                    # 场高斯求积点
                    rgi =   (gi != gj ? getGQPTriSglr(tri, gi) : rgj)
                    # ρt
                    ρt  =   rgi - freeVm
                    # 计算结果，利用对称性加速计算，因此要对非对角线结果乘2
                    Ztemp   +=  (gi != gj ? 2*(ρt ⋅ ρs - C4divk²)*gstarw[gi, gj] : (ρt ⋅ ρs - C4divk²)*gstarw[gi, gj])
                end #for gi
            end #for gj
            
            # 奇异部分
            if  mi == ni
                # 找出三边
                a, b, c     =   abs.(tri.edgel[MoM_Basics.Vec3IdxCircle[mi:mi+2]])
                #奇异项计算
                Ztemp     +=   singularF21(a, b, c, tri.area^2) - F1
            else
                # 找出三边
                local temp  =   6 - mi - ni
                a, b, c     =   abs.(tri.edgel[MoM_Basics.Vec3IdxCircle[temp:temp+2]])
                # 奇异项计算
                Ztemp     +=   singularF22(a, b, c, tri.area^2) - F1
            end #if


            Ztemp   *=  lm*ln*JKηdiv16π
            # 将结果写入目标数组
            Ztri[mi, ni]    =   Ztemp
            Ztri[ni, mi]    =   Ztemp

        end #for mi

    end # for ni 

    return Ztri

end

"""
本函数用于计算PEC的EFIE阻抗矩阵。
输入信息：
trianglesInfo:  为包含三角形信息实例的向量
nrwg        :  基函数数目
返回值
Zmat         :  阻抗矩阵

注意，此程序由于采用的对三角形循环计算，因此在并行化时，会出现不同线程计算出同一个矩阵元，导致写入冲突，因此要加线程锁保证结果写入正确
"""
function impedancemat4EFIE4PEC(trianglesInfo::Vector{TriangleInfo{IT, FT}}, nrwg::Integer, bfT::Type{BFT}) where {IT, FT, BFT<:RWG}
    # 初始化阻抗矩阵
    Zmat    =   zeros(Complex{FT}, (nrwg, nrwg))
    impedancemat4EFIE4PEC!(Zmat, trianglesInfo, bfT)
    return Zmat
end


"""
本函数用于在有矩阵的情况下计算PEC的EFIE阻抗矩阵。
输入信息：
Zmat
trianglesInfo:  为包含三角形信息实例的向量
nrwg        :  基函数数目
返回值
Zmat         :  阻抗矩阵

注意，此程序由于采用的对三角形循环计算，因此在并行化时，会出现不同线程计算出同一个矩阵元，导致写入冲突，因此要加线程锁保证结果写入正确

"""
function impedancemat4EFIE4PEC!(Zmat::Matrix{Complex{FT}}, trianglesInfo::Vector{TriangleInfo{IT, FT}}, ::Type{BFT}) where {IT, FT, BFT<:RWG}
    
    # 索引
    trisIdx =   eachindex(trianglesInfo)
    # 常数
    Rsglr       =   Params.Rsglr
    # 三角形数
    trisnum =   trisIdx.stop
    # 线程锁防止对同一数据写入出错
    lockZ   =   SpinLock()
    # 矩阵大小
    nbf     =   size(Zmat, 1)
    # Progress Meter
    pmeter  =   Progress(trisnum; desc = "Calculating Z (RWG, EFIE) ($nbf × $nbf)")

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
                Ztri  =  EFIEOnTris(trit)
                # 写入数据
                lock(lockZ)
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = trit.inBfsID[ni]
                    # 判断边是不是基函数（边缘不算）
                    (m == 0 || n == 0) && continue
                    # 往矩阵填充结果
                    Zmat[m, n] += Ztri[mi, ni]
                end
                unlock(lockZ)

            elseif Rts < Rsglr
                # 需要进行近奇异性处理的场源三角形
                Ztri    =   EFIEOnNearTris(trit, tris)
                
                # 写入数据
                lock(lockZ)
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = tris.inBfsID[ni]

                    # 判断边是不是基函数（边缘不算）
                    (m == 0 || n == 0) && continue
                    
                    Zmat[m, n] += Ztri[mi, ni]
                    Zmat[n, m] += Ztri[mi, ni]
                end
                unlock(lockZ)

            else
                # 正常高斯求积
                # 计算三角形相关的(3*3)个矩阵元的结果
                Ztri    =   EFIEOnTris(trit, tris)
                
                # 写入数据
                lock(lockZ)
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m = trit.inBfsID[mi]
                    n = tris.inBfsID[ni]

                    # 判断边是不是基函数（边缘不算）
                    (m == 0 || n == 0) && continue
                    
                    Zmat[m, n] += Ztri[mi, ni]
                    Zmat[n, m] += Ztri[mi, ni]
                end
                unlock(lockZ)

            end # if

        end #for trisj

        next!(pmeter)

    end #for triti

    return Zmat
    
end