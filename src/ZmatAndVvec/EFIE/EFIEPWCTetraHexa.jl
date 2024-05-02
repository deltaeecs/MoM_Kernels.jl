
"""
采用 PWC 基函数
计算六面体和四面体上相关的9个阻抗矩阵元，
此函数方法用于计算场源六面体不重合且相隔较远的情况，因此输入有两个六面体信息类型实例
输入：
hexat::HexahedraInfo,   场六面体
tetras::TetrahedraInfo  源四面体
计算：
`jk₀η₀∫ₜ∫ₛ(I + 1/k²∇∇)G(R)dV′dV`
注意为方便对称性快速填充矩阵元，没有加入 κ 项，因此后续填充时要注意加上。
"""
function EFIEOnHexaTetraPWC(hexat::HexahedraInfo{IT, FT, CT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    
    # 场源求积点
    rgt     =   getGQPHexa(hexat)
    rgs     =   getGQPTetra(tetras)

    # 采样点对应的体元
    dVtdVs  =   hexat.volume*tetras.volume
    # 常数项
    JK_0    =   Params.JK_0
    Jη_0divKdVtdVs  =   Params.Jη_0divK * dVtdVs
    k²      =   Params.k²
    # 距离向量
    Rtsvec  =   zero(MVec3D{FT})

    # 储存结果的临时数组
    # re = zero(MMatrix{3, 3, CT})
    ## 计算矩阵元并矢
    # 对源求积点循环
    @inbounds for gj in 1:GQPNTetra
        # 源高斯求积点
        rgj  =   view(rgs, :, gj)
        # 对场求积点循环
        for gi in 1:GQPNHexa
            # 场高斯求积点
            rgi  =  view(rgt, :, gi)

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
            GR  =   exp(-JK_0*Rts)*div4π*divR*HexaGQInfo.weight[gi]*TetraGQInfo.weight[gj]
            # 计算矩阵元并叠加
            Zts .+=     (Jη_0divKdVtdVs * GR) * ((I - R̂R̂) * k² - (I/3 -  R̂R̂) * (3*jkplusR1stdivR1st) )
            
        end # for gi
    end #for gj

    return Zts
end

"""
采用 PWC 基函数
计算六面体和四面体上相关的9个阻抗矩阵元，
此函数方法用于计算场源六面体不重合且相隔较远的情况，因此输入有两个六面体信息类型实例
输入：
tetrat::TetrahedraInfo  场四面体
hexas::HexahedraInfo,   源六面体
计算：
`jk₀η₀∫ₜ∫ₛ(I + 1/k²∇∇)G(R)dV′dV`
注意为方便对称性快速填充矩阵元，没有加入 κ 项，因此后续填充时要注意加上。
"""
function EFIEOnHexaTetraPWC(tetrat::TetrahedraInfo{IT, FT, CT}, hexas::HexahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    Zst     =   EFIEOnHexaTetraPWC(hexas, tetrat)
    Zts    .=   transpose(Zst)

    return Zts
end


"""
采用 PWC 基函数
计算六面体上相关的9个阻抗矩阵元，
此函数方法用于计算场源六面体不重合且相隔较近的情况，因此输入有两个六面体信息类型实例
输入：
hexat  hexas     :   HexahedraInfo, 场六面体和六面体
计算：
jkη₀∫ₜ∫ₛ(I + 1/k²∇∇)G(R)dV′dV
其中， ∫ₜ∫ₛ∇∇G(R)dV′dV = ∫ₜ∑ᵢn̂ᵢ(∫ᵢR̂(jk + 1/R)G(R)dS′)dV
计算得到结果为并矢::
jη₀/k ∫∫ (k²I + ∇∇)G(R) dV'dV
Kᵣⁿ  =   ∫ Rⁿ dV'
K̂ᵣⁿ  =   ∫ R̂Rⁿ dV'
注意为方便对称性快速填充矩阵元，没有加入 κ 项，因此后续填充时要注意加上。
"""
function EFIEOnNearHexaTetraPWC(hexat::HexahedraInfo{IT, FT, CT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    
    # 场求积点
    rgt     =   getGQPHexaSglr(hexat)

    # 采样点对应的体元
    dVt  =   hexat.volume
    # 常数项
    Jη_0divKdVt =   Params.Jη_0divK * dVt
    ## 计算矩阵元并矢
    # 对场求积点循环
    @inbounds for gi in 1:GQPNHexaSglr
        # 源高斯求积点
        rgi  =      rgt[:, gi]
        # 计算 L 算子并矢并乘以权重
        Zts .+=     HexaGQInfoSglr.weight[gi] .* volumeSingularityLOpDyad(rgi, tetras)
    end

    # 补上常数项
    Zts .*= Jη_0divKdVt

    return Zts
end

"""
采用 PWC 基函数
计算六面体和四面体上相关的9个阻抗矩阵元，
此函数方法用于计算场源六面体不重合且相隔较远的情况，因此输入有两个六面体信息类型实例
输入：
tetrat::TetrahedraInfo  场四面体
hexas::HexahedraInfo,   源六面体
计算：
`jk₀η₀∫ₜ∫ₛ(I + 1/k²∇∇)G(R)dV′dV`
注意为方便对称性快速填充矩阵元，没有加入 κ 项，因此后续填充时要注意加上。
"""
function EFIEOnNearHexaTetraPWC(tetrat::TetrahedraInfo{IT, FT, CT}, hexas::HexahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    Zst     =   EFIEOnNearHexaTetraPWC(hexas, tetrat)
    Zts    .=   transpose(Zst)

    return Zts
end

"""
本函数用于计算介质体的 PWC 基函数下的 EFIE 阻抗矩阵。
输入信息：
hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}},  为包含六面体信息实例的向量
tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, 为包含四面体信息实例的向量
nPWC        :  基函数数目
返回值
Zmat        :  阻抗矩阵
"""
function impedancemat4VIE!(Zmat::Matrix{CT}, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}},
                            tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, ::Type{BFT};
                            discreteVar = SimulationParams.discreteVar) where {IT, FT, CT, BFT<:PWC}
    
    # 六面体数
    hexasnum    =   length(hexasInfo)
    # 判断体电流的离散方式
    discreteJ::Bool = discreteVar == "J"
    # 常数
    Rsglr   =   Params.Rsglr
    # Progress Meter
    nbf     =   size(Zmat, 1)
    pmeter  =   Progress(hexasnum; desc = "Calculating Z (PWC)($nbf × $nbf)...")
    # 外层定义为场基函数循环
    @threads for hexat in hexasInfo
        # 场六面体介质对比度
        κₜ  =   hexat.κ

        Rsglrlc =   Rsglr/sqrt(norm(hexat.ε)/ε_0)
        @inbounds for tetras in tetrasInfo
            # 源六面体介质对比度
            κₛ  =   tetras.κ
            # 场源距离
            local Rts   =   dist(hexat.center, tetras.center)
            # isapprox(Rts, Rsglrlc, rtol = 1e-2) && @show ti, sj
            # 判断二者远近，调用不同精度的矩阵元处理函数
            if Rts < Rsglr
                # 需要进行近奇异性处理的场源六面体
                Zts    =   EFIEOnNearHexaTetraPWC(hexat, tetras)
                # 写入数据，利用对称性填充，因此要避免重合时重复填充
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m   =   hexat.inBfsID[mi]
                    n   =   tetras.inBfsID[ni]
                    # 写入
                    if discreteJ
                        Zmat[m, n]  =   Zts[mi, ni]
                        Zmat[n, m]  =   Zts[mi, ni]
                    else
                        Zmat[m, n]  =   Zts[mi, ni]*κₛ
                        Zmat[n, m]  =   Zts[mi, ni]*κₜ
                    end
                end
            else
                # 正常高斯求积
                # 计算六面体相关的(3*3)个矩阵元的结果
                Zts    =   EFIEOnHexaTetraPWC(hexat, tetras)
                
                # 写入数据
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m   =   hexat.inBfsID[mi]
                    n   =   tetras.inBfsID[ni]
                    # 写入
                    if discreteJ
                        Zmat[m, n]  =   Zts[mi, ni]
                        Zmat[n, m]  =   Zts[mi, ni]
                    else
                        Zmat[m, n]  =   Zts[mi, ni]*κₛ
                        Zmat[n, m]  =   Zts[mi, ni]*κₜ
                    end
                end

            end # if

        end #for sj

        next!(pmeter)

    end #for ti

    return nothing
    
end

"""
本函数用于计算介质体的 PWC 基函数下的 EFIE 阻抗矩阵。
输入信息：
hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}},  为包含六面体信息实例的向量
tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, 为包含四面体信息实例的向量
nPWC        :  基函数数目
返回值
Zmat        :  阻抗矩阵
"""
function impedancemat4VIE!(Zmat::Matrix{CT}, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}},
                            hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}, bfT::Type{BFT}) where {IT, FT, CT, BFT<:PWC}
    impedancemat4VIE!(Zmat, hexasInfo, tetrasInfo, bfT)

end

"""
本函数用于计算混合网格（四面体+六面体）下介质体的 PWC 基函数下的 EFIE 阻抗矩阵。
输入信息：
geosInfo    :  为包含四面体信息实例的向量
nPWC        :  基函数数目
返回值
Zmat        :  阻抗矩阵
"""
function impedancemat4VIE(geosInfoV::AbstractVector{VT}, nPWC::Integer, bfT::Type{BFT}) where {VT<:AbstractVector, BFT<:PWC}
    CT = typeof(geosInfoV[1][1].ε)
    # 初始化阻抗矩阵
    Zmat    =   zeros(CT, (nPWC, nPWC))
    # 计算各自阻抗矩阵
    for geosInfo in geosInfoV
        impedancemat4VIE!(Zmat, geosInfo, bfT)
    end
    # 计算相互之间的阻抗矩阵
    for i in 1:length(geosInfoV)
        geosA = geosInfoV[i]
        for j in (i+1):length(geosInfoV)
            geosB = geosInfoV[j]
            impedancemat4VIE!(Zmat, geosA, geosB, bfT)
        end
    end

    return Zmat
    
end