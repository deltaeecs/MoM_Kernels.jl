
"""
计算三角形上相关9个阻抗矩阵元，
此函数方法用于计算场源三角形不重合且相隔较远的情况，因此输入有两个个三角形信息类型实例
输入
trit， tris     :   TriangleInfo, 场三角形和源三角形
"""
function MFIEOnTris(trit::TriangleInfo{IT, FT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
    CT  =   Complex{FT}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    Zst     =   zeros(CT, 3, 3)
    # 常数
    JK_0    =   Params.JK_0
    n̂t      =   trit.facen̂
    n̂s      =   tris.facen̂
    # 预分配内存
    gw      =   zero(MMatrix{GQPNTri, GQPNTri, CT})
    # gw      =   zeros(CT, GQPNTri, GQPNTri)
    
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

    # 预分配 ρs ρt
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
            let Zmn = zero(CT), Znm = zero(CT)
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
                        # 计算 Zmn\Znm
                        temp    =   (JK_0 + divr) * divr * gw[gi, gj]
                        # Zmn    +=   (ρmi × n̂t) ⋅ (rvec × ρnj) * temp
                        # Znm    -=   (ρnj × n̂s) ⋅ (rvec × ρmi) * temp
                        ρmiρnj  =   ρmi ⋅ ρnj
                        Zmn    +=   ((ρmi ⋅ rvec) * (n̂t ⋅ ρnj) - (n̂t ⋅ rvec) * ρmiρnj) * temp
                        Znm    -=   ((ρnj ⋅ rvec) * (n̂s ⋅ ρmi) - (n̂s ⋅ rvec) * ρmiρnj) * temp
                    end #for gi
                end #for gj
                temp    =   lm*ln*ηdiv16π
                Zmn    *=   temp
                Znm    *=   temp
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
function MFIEOnNearTris(trit::TriangleInfo{IT, FT}, tris::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
    CT  =   Complex{FT}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    Zst     =   zeros(CT, 3, 3)
    # 常数
    JK_0    =   Params.JK_0
    n̂t      =   trit.facen̂
    n̂s      =   tris.facen̂
    # 预分配内存
    gw      =   zero(MMatrix{GQPNTriSglr, GQPNTriSglr, CT})
    # gw      =   zeros(CT, GQPNTri, GQPNTri)
    
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

    # 预分配 ρs ρt
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
            let Zmn = zero(CT), Znm = zero(CT)
                # 非奇异直接高斯求积
                for gj in 1:GQPNTriSglr
                    # 源高斯求积点
                    rgj     =   view(rgs, :, gj)
                    # ρn
                    ρnj    .=   rgj .- freeVn
                    # 对场求积点循环
                    for gi in 1:GQPNTriSglr
                        # 场高斯求积点
                        rgi     =   view(rgt, :, gi)
                        # ρm
                        ρmi    .=   rgi .- freeVm
                        # rvec
                        rvec   .=   rgi .- rgj
                        # r
                        divr    =   1/norm(rvec)
                        # 计算 Zmn\Znm
                        Zmn    +=   (ρmi × n̂t) ⋅ (rvec × ρnj) * ((JK_0 + divr) * divr * gw[gi, gj])
                        Znm    -=   (ρnj × n̂s) ⋅ (rvec × ρmi) * ((JK_0 + divr) * divr * gw[gi, gj])                  
                    end #for gi
                end #for gj
                temp    =   lm*ln*ηdiv16π
                Zmn    *=   temp
                Znm    *=   temp
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
此函数方法用于计算场源三角形重合的情况，因此输入只有一个三角形信息类型实例
输入
tri     :   TriangleInfo
"""
function MFIEOnTris(tri::TriangleInfo{IT, FT}) where {IT<: Integer, FT<:AbstractFloat}
    CT  =   Complex{FT}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)

    # 高斯求积点
    rgt     =   getGQPTri(tri)

    # 常数
    ηdiv8S  =   η_0/(8tri.area)
    
    # 预分配 ρs ρt
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
            Ztemp   =   zero(CT)
            # 非奇异部分
            for gj in 1:GQPNTri
                # 源高斯求积点
                rgj =   view(rgt, :, gj)
                # ρs
                ρnj .=   rgj .- freeVn
                ρmi .=   rgj .- freeVm
                # 计算
                Ztemp +=    ρnj ⋅ ρmi * TriGQInfo.weight[gj]

            end #for gj

            Ztemp   *=  lm*ln*ηdiv8S
            # 将结果写入目标数组
            Zts[mi, ni] =   Ztemp
            Zts[ni, mi] =   Ztemp

        end #for mi
    end # for ni 

    return Zts

end

"""本函数用于计算PEC的MFIE阻抗矩阵。
输入信息：
trianglesInfo:  为包含三角形信息实例的向量
nrwg        :  基函数数目
返回值
Zmat         :  阻抗矩阵

注意，此程序由于采用的对三角形循环计算，因此在并行化时，会出现不同线程计算出同一个矩阵元，导致写入冲突，因此要加线程锁保证结果写入正确

"""
function impedancemat4MFIE4PEC(trianglesInfo::Vector{TriangleInfo{IT, FT}}, nrwg::Integer, ::Type{BFT}) where {IT, FT, BFT<:RWG}
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
                Zts  =  MFIEOnTris(trit)
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
                Zts, Zst    =   MFIEOnNearTris(trit, tris)
                
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
                Zts, Zst    =   MFIEOnTris(trit, tris)
                
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
