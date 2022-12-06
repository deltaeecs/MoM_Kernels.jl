"""
计算平面波在 RWG 基函数上的激励向量
输入：
source          ::ST, 波源
trianglesInfo   ::Vector{TriangleInfo{IT, FT}}，保存三角形信息的向量
nbf             ::Integer，基函数数目
"""
function excitationVectorCFIE!(V::Vector{Complex{FT}}, source::ST, trianglesInfo::Vector{TriangleInfo{IT, FT}}, sbfType = VSBFTypes.sbfType) where{ST<:ExcitingSources, IT<:Integer, FT<:Real}
    
    nt  =   length(trianglesInfo)
    pmeter  =   Progress(nt, "Calculating V...")
    lockV   =   SpinLock()
    # 开始对四面体形循环计算
    @threads for i in 1:nt
        next!(pmeter)
        tri = trianglesInfo[i]
        # 三角形相关激励向量
        Vtri    =   excitationVectorCFIE(source, tri, sbfType)
        # 写入结果
        for ni in 1:3
            n = tri.inBfsID[ni]
            # 跳过半基函数
            (n == 0) && continue
            lock(lockV)
            V[n] += Vtri[ni]
            unlock(lockV)
        end
    end # for tri

    nothing
end #for function


"""
计算平面波在 RWG 基函数上的激励向量
输入：
source          ::ST, 波源
trianglesInfo   ::Vector{TriangleInfo{IT, FT}}，保存三角形信息的向量
nbf             ::Integer，基函数数目
"""
function excitationVectorCFIE(source::ST, trianglesInfo::Vector{TriangleInfo{IT, FT}}, nbf::Integer) where{ST<:ExcitingSources, IT<:Integer, FT<:Real}

    # 预构造激励向量
    V   =   zeros(Complex{FT}, nbf)
    # 面基函数的类型
    sbfType =   VSBFTypes.sbfType
    excitationVectorCFIE!(V, source, trianglesInfo, sbfType)
    # 返回
    return V
end #for function

"""
计算平面波在给定三角形的三个 半RWG 基函数上的激励向量
输入：
source  ::ST, 波源
tri     ::TriangleInfo{IT, FT}，三角形信息
"""
function excitationVectorCFIE(source::ST, tri::TriangleInfo{IT, FT}, ::Type{BFT}) where {ST<:ExcitingSources, IT<:Integer, FT<:Real, BFT<:RWG}
    # 复数类型
    CT  =   Complex{FT}
    # 混合场积分方程系数
    α   =   Params.CFIEα
    # 所有的高斯求积点
    rvecgs  =   getGQPTri(tri)
    # 高斯求积点上的电场值
    Egs     =   zero(MMatrix{3, GQPNTri, CT})
    for  g in 1:GQPNTri
        Egs[:, g]    =  sourceEfield(source, rvecgs[:, g])
    end    
    # 高斯求积点上的磁场值
    Hgs     =   zero(MMatrix{3, GQPNTri, CT})
    for  g in 1:GQPNTri
        Hgs[:, g]    =  sourceHfield(source, rvecgs[:, g])
    end
    # 面的外法向量
    n̂   =   tri.facen̂
    # 预定义并置零结果
    Vtri =   zero(MVec3D{CT})
    # 对三角形的三个基函数循环
    @views for mi in 1:3

        # 所在基函数id
        m   =   tri.inBfsID[mi]
        # 判断边是不是基函数（边缘不算）
        m == 0 && continue
        # 对高斯求积点循环
        VE  =   zero(CT)
        VM  =   zero(CT)
        for g in 1:GQPNTri           
            # ρm
            ρm  =   rvecgs[:, g] - tri.vertices[:, mi]
            # 累加结果
            VM +=   (ρm ⋅ (n̂ × view(Hgs, :, g)))*TriGQInfo.weight[g]
            VE +=   (ρm ⋅ view(Egs, :, g))*TriGQInfo.weight[g]
        end
        # 最后的计算
        Vtri[mi]    =   α*VE + (1-α)*η_0*VM
        Vtri[mi]   *=  tri.edgel[mi]/2
    end

    return Vtri

end
