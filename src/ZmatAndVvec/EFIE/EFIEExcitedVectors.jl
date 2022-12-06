"""
计算平面波在 RWG 基函数上的激励向量
输入：
source          ::ST, 波源
trianglesInfo   ::Vector{TriangleInfo{IT, FT}}，保存三角形信息的向量
nbf             ::Integer，基函数数目
"""
function excitationVectorEFIE!(V::Vector{Complex{FT}}, source::ST, trianglesInfo::Vector{TriangleInfo{IT, FT}}, sbfType = VSBFTypes.sbfType) where{ST<:ExcitingSources, IT<:Integer, FT<:Real}
    
    nt  =   length(trianglesInfo)
    pmeter  =   Progress(nt, "Calculating V...")
    lockV   =   SpinLock()
    # 开始对四面体形循环计算
    @threads for i in 1:nt
        next!(pmeter)
        tri = trianglesInfo[i]
        # 三角形相关激励向量
        Vtri    =   excitationVectorEFIE(source, tri, sbfType)
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
function excitationVectorEFIE(source::ST, trianglesInfo::Vector{TriangleInfo{IT, FT}}, nbf::Integer) where{ST<:ExcitingSources, IT<:Integer, FT<:Real}

    # 预构造激励向量
    V   =   zeros(Complex{FT}, nbf)
    # 面基函数的类型
    sbfType =   VSBFTypes.sbfType
    excitationVectorEFIE!(V, source, trianglesInfo, sbfType)

    # 返回
    return V
end #for function



"""
计算平面波在给定三角形的三个 半RWG 基函数上的激励向量
输入：
source  ::ST, 波源
tri     ::TriangleInfo{IT, FT}，三角形信息
"""
function excitationVectorEFIE(source::ST, tri::TriangleInfo{IT, FT}, ::Type{BFT}) where {ST<:ExcitingSources, IT<:Integer, FT<:Real, BFT<:RWG}
    # 复数类型
    CT  =   Complex{FT}
    # 所有的高斯求积点
    rvecgs  =   getGQPTri(tri)
    # 高斯求积点上的场值
    Egs     =   zero(MMatrix{3, GQPNTri, CT})
    for  g in 1:GQPNTri
        Egs[:, g]    =  sourceEfield(source, rvecgs[:, g])
    end
    # 预定义并置零结果
    Vtri =   zero(MVec3D{CT})
    # 对三角形的三个基函数循环
    @views for mi in 1:3

        # 所在基函数id
        m       =   tri.inBfsID[mi]
        # 判断边是不是基函数（边缘不算）
        m == 0 && continue
        # 对高斯求积点循环
        for g in 1:GQPNTri           
            # ρm
            ρm      =   rvecgs[:, g] - tri.vertices[:, mi]
            # 累加结果
            Vtri[mi]  +=   (ρm ⋅ view(Egs, :, g))*TriGQInfo.weight[g]
        
        end

        # 最后的计算
        Vtri[mi]   *=  tri.edgel[mi]/2
    end

    return Vtri

end

"""
计算平面波在 基函数 上的激励向量
输入：
source          ::ST, 平面波源
tetrasInfo      ::Vector{TetrahedraInfo{IT, FT, CT}}，保存四面体信息的向量
nbf             ::Integer，基函数数量  
"""
function excitationVectorEFIE(source::ST, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, nbf::Integer, vbfType =   VSBFTypes.vbfType) where {ST<:ExcitingSources, IT<:Integer, FT<:Real, CT<:Complex{FT}}

    # 预构造激励向量
    V   =   zeros(Complex{FT}, nbf)
    # 开始对四面体形循环计算
    excitationVectorEFIE!(V, source, tetrasInfo, vbfType)

    # 返回
    return V
end # function

"""
计算平面波在 基函数 上的激励向量
输入：
source          ::ST, 平面波源
tetrasInfo      ::Vector{TetrahedraInfo{IT, FT, CT}}，保存四面体信息的向量
nbf             ::Integer，基函数数量  
"""
function excitationVectorEFIE!(V::Vector{Complex{FT}}, source::ST, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, vbfType =   VSBFTypes.vbfType) where {ST<:ExcitingSources, IT<:Integer, FT<:Real, CT<:Complex{FT}}

    # print("Calculating Excitation Vector...")
    nt  =   length(tetrasInfo)
    pmeter  =   Progress(nt, "Calculating V...")
    lockV   =   SpinLock()
    # 开始对四面体形循环计算
    @threads for i in eachindex(tetrasInfo)
        next!(pmeter)
        tetra = tetrasInfo[i]
        # 四面体相关激励向量
        Vtetra    =   excitationVectorEFIE(source, tetra, vbfType)
        # 写入结果
        lock(lockV)
        V[tetra.inBfsID] .+= Vtetra
        unlock(lockV)
    end # for tetra

    # 返回
    return V
end # function

"""
计算平面波在给定四面体的四个 半SWG 基函数上的激励向量
输入：
source      ::ST, 波源
tetra       ::TetrahedraInfo{IT, FT, CT}，四面体信息
"""
function excitationVectorEFIE(source::ST, tetra::TetrahedraInfo{IT, FT, CT}, ::Type{BFT}) where {ST<:ExcitingSources, IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:LinearBasisFunction}
    # 所有的高斯求积点
    rvecgs  =   getGQPTetra(tetra)
    # 高斯求积点上的场值
    Egs     =   zero(MMatrix{3, GQPNTetra, CT})
    for ig in 1:GQPNTetra
        Egs[:, ig]    =  sourceEfield(source, rvecgs[:, ig])
    end
    # 预定义并置零结果
    Vtetra =   zero(MVector{4, CT})
    # 对四面体的四个 SWG 基函数循环
    @views for mi in 1:4
        # 对高斯求积点循环
        for g in 1:GQPNTetra           
            # ρm
            ρm      =   rvecgs[:, g] .- tetra.vertices[:, mi]
            # 累加结果
            Vtetra[mi]  +=   (ρm ⋅ view(Egs, :, g))*TetraGQInfo.weight[g]
        end
        # 最后的计算
        Vtetra[mi]   *=  tetra.facesArea[mi]/3
    end

    return Vtetra

end

"""
计算平面波在给定四面体的三个 PWC 基函数上的激励向量
输入：
source      ::ST, 波源
tetra       ::TetrahedraInfo{IT, FT, CT}，四面体信息
"""
function excitationVectorEFIE(source::ST, tetra::TetrahedraInfo{IT, FT, CT}, ::Type{BFT}) where {ST<:ExcitingSources, IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:ConstBasisFunction}
    # 所有的高斯求积点
    rvecgs  =   getGQPTetra(tetra)
    # 预定义并置零结果
    Vtetra =   zero(MVec3D{CT})
    # 对高斯求积点循环
    for ig in 1:GQPNTetra
        # 高斯求积点上的场值
        Eg  =  sourceEfield(source, rvecgs[:, ig])
        # 累加结果
        Vtetra  .+=   Eg * TetraGQInfo.weight[ig]
    end
    # 最后的计算
    Vtetra   .*=  tetra.volume

    return Vtetra

end


"""
计算平面波在 基函数 上的激励向量
输入：
source          ::ST, 平面波源
hexasInfo      ::Vector{HexahedraInfo{IT, FT, CT}}，保存六面体信息的向量
nbf             ::Integer，基函数数量  
"""
function excitationVectorEFIE(source::ST, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}, nbf::Integer, vbfType =   VSBFTypes.vbfType) where{ST<:ExcitingSources, IT<:Integer, FT<:Real, CT<:Complex{FT}}

    # print("Calculating Excitation Vector...")
    # 预构造激励向量
    V   =   zeros(Complex{FT}, nbf)
    # 开始对六面体形循环计算
    excitationVectorEFIE!(V, source, hexasInfo, vbfType)

    # 返回
    return V
end # function

"""
计算平面波在 基函数 上的激励向量
输入：
source          ::ST, 平面波源
hexasInfo       ::Vector{HexahedraInfo{IT, FT, CT}}，保存六面体信息的向量
nbf             ::Integer，基函数数量  
"""
function excitationVectorEFIE!(V::Vector{CT}, source::ST, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}, vbfType =   VSBFTypes.vbfType) where{ST<:ExcitingSources, IT<:Integer, FT<:Real, CT<:Complex{FT}}

    # print("Calculating Excitation Vector...")
    n   =   length(hexasInfo)
    pmeter  =   Progress(n, "Calculating V...")
    lockV   =   SpinLock()
    # 开始对四面体形循环计算
    @threads for i in eachindex(hexasInfo)
        next!(pmeter)
        hexa = hexasInfo[i]
        # 六面体相关激励向量
        Vhexa    =   excitationVectorEFIE(source, hexa, vbfType)
        # 写入结果
        lock(lockV)
        V[hexa.inBfsID] .+= Vhexa
        unlock(lockV)
    end # for tetra

    # 返回
    return V
end # function


"""
计算平面波在给定六面体的六个 半SWG 基函数上的激励向量
输入：
source      ::ST, 波源
hexa        ::HexahedraInfo{IT, FT, CT}，六面体信息
"""
function excitationVectorEFIE(source::ST, hexa::HexahedraInfo{IT, FT, CT}, ::Type{BFT}) where {ST<:ExcitingSources, IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:LinearBasisFunction}
    # 所有的高斯求积点
    rvecgs  =   getGQPHexa(hexa)
    # 高斯求积点上的场值
    Egs     =   zero(MMatrix{3, GQPNHexa, CT})
    for g in 1:GQPNHexa
        Egs[:, g]    =  sourceEfield(source, rvecgs[:, g])
    end
    # 预定义并置零结果
    Vhexa =   zero(MVector{6, CT})
    # 对六面体六个基函数循环计算
    @views for mi in 1:6
        freeVms     =   getFreeVns(hexa, mi)
        # 对高斯求积点循环
        for g in 1:GQPNHexa
            # 计算场六面体的 “自由端” id
            idm3D   =   GQ1DID2GQ3DIDVector[g]
            # 计算源六面体的 “自由端” id
            idm     =   getFreeVIDFromGQ3DID(idm3D, mi)
            # ρm
            ρm      =   rvecgs[:, g] .- view(freeVms, :, idm)
            # 累加结果
            Vhexa[mi]  +=   (ρm ⋅ view(Egs, :, g))*HexaGQInfo.weight[g]
        end
        # 最后的计算
        Vhexa[mi]   *=  hexa.facesArea[mi]
    end

    return Vhexa

end

"""
计算平面波在给定四面体的三个 PWC 基函数上的激励向量
输入：
source      ::ST, 波源
hexa        ::HexahedraInfo{IT, FT, CT}，六面体信息
"""
function excitationVectorEFIE(source::ST, hexa::HexahedraInfo{IT, FT, CT}, ::Type{BFT}) where {ST<:ExcitingSources, IT<:Integer, FT<:Real, CT<:Complex{FT}, BFT<:ConstBasisFunction}
    # 所有的高斯求积点
    rvecgs  =   getGQPHexa(hexa)
    # 预定义并置零结果
    Vhexa =   zero(MVec3D{CT})
    # 对高斯求积点循环
    for ig in 1:GQPNHexa
        # 高斯求积点上的场值
        Eg  =  sourceEfield(source, rvecgs[:, ig])
        # 累加结果
        Vhexa  .+=   Eg * HexaGQInfo.weight[ig]
    end
    # 最后的计算
    Vhexa   .*=  hexa.volume

    return Vhexa

end

using MoM_Basics:getBFTfromCellT
"""
计算平面波在 基函数 上的激励向量
输入：
source          ::ST, 平面波源
tetrasInfo      ::Vector{TetrahedraInfo{IT, FT, CT}}，保存六面体信息的向量
nbf             ::Integer，基函数数量  
"""
function excitationVectorEFIE(source::ST, geosInfo::Vector{VT}, nbf::Integer) where{ST<:ExcitingSources, VT<:AbstractVector}

    # 预构造激励向量
    V   =   zeros(typeof(first(first(geosInfo)).ε), nbf)
    # 开始对网格类型（三角形、四边形）循环计算
    for geo in geosInfo
        # 网格类型相关激励向量
        excitationVectorEFIE!(V, source, geo, getBFTfromCellT(eltype(geo)))
    end # for geo
    # 返回
    return V
end # function