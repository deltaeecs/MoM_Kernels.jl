
"""
采用 PWC 基函数
计算四面体上相关的9个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较远的情况，因此输入有两个四面体信息类型实例
输入：
tetrat  tetras     :   TetrahedraInfo, 场四面体和四面体
计算：
jk₀η₀∫ₜ∫ₛ(I + 1/k²∇∇)G(R)dV′dV
注意为方便对称性快速填充矩阵元，没有加入 κ 项，因此后续填充时要注意加上。
"""
function EFIEOnTetrasPWC(tetrat::TetrahedraInfo{IT, FT, CT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    
    # 场源求积点
    rgt     =   getGQPTetra(tetrat)
    rgs     =   getGQPTetra(tetras)

    # 采样点对应的体元
    dVtdVs  =   tetrat.volume*tetras.volume
    # 常数项
    JK_0    =   Params.JK_0
    k²      =   Params.k² 
    Jη_0divKdVtdVs  =   Params.Jη_0divK * dVtdVs
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
        for gi in 1:GQPNTetra
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
            GR  =   exp(-JK_0*Rts)*div4π*divR*TetraGQInfo.weight[gi]*TetraGQInfo.weight[gj]

            # 计算矩阵元并叠加
            # Zts .+=     (Jη_0divKdVtdVs * GR) * ((I - R̂R̂) * k² - (I/3 -  R̂R̂) * (3*jkplusR1stdivR1st) )
            for j in 1:3, i in 1:3
                i == j ? begin
                    Zts[i, j]   +=   (Jη_0divKdVtdVs * GR) * ((1 - R̂R̂[i, j]) * k² - (1 - 3R̂R̂[i, j])*jkplusR1stdivR1st )
                end : begin
                    Zts[i, j]   +=   (Jη_0divKdVtdVs * GR) * ( - R̂R̂[i, j] * k² + 3R̂R̂[i, j]*jkplusR1stdivR1st )
                end
            end

        end # for gi
    end #for gj

    return Zts
end

"""
采用 PWC 基函数
计算四面体上相关的9个阻抗矩阵元，
此函数方法用于计算场源四面体不重合且相隔较近的情况，因此输入有两个四面体信息类型实例
输入：
tetrat  tetras     :   TetrahedraInfo, 场四面体和四面体
计算：
jkη₀∫ₜ∫ₛ(I + 1/k²∇∇)G(R)dV′dV
其中， ∫ₜ∫ₛ∇∇G(R)dV′dV = ∫ₜ∑ᵢn̂ᵢ(∫ᵢR̂(jk + 1/R)G(R)dS′)dV
计算得到结果为并矢::
jη₀/k ∫∫ (k²I + ∇∇)G(R) dV'dV
Kᵣⁿ  =   ∫ Rⁿ dV'
K̂ᵣⁿ  =   ∫ R̂Rⁿ dV'
注意为方便对称性快速填充矩阵元，没有加入 κ 项，因此后续填充时要注意加上。
"""
function EFIEOnNearTetrasPWC(tetrat::TetrahedraInfo{IT, FT, CT}, tetras::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    
    # 场求积点
    rgt     =   getGQPTetraSglr(tetrat)

    # 采样点对应的体元
    dVt  =   tetrat.volume
    # 常数项
    Jη_0divKdVt =   Params.Jη_0divK * dVt
    ## 计算矩阵元并矢
    # 对场求积点循环
    @inbounds for gi in 1:GQPNTetraSglr
        # 场高斯求积点
        rgi     =   view(rgt, :, gi)
        # 计算 L 算子并矢并乘以权重
        Zts   .+=   TetraGQInfoSglr.weight[gi] .* volumeSingularityLOpDyad(rgi, tetras)
    end

    # 补上常数项
    Zts .*= Jη_0divKdVt

    return Zts
end

"""
采用 PWC 基函数
计算四面体上相关的9个阻抗矩阵元，
此函数方法用于计算场源四面体重合的情况，因此输入有一个四面体信息类型实例
输入：
tetrat   TetrahedraInfo, 场四面体和四面体
计算：
jkη₀∫ₜ∫ₛ(I + 1/k²∇∇)G(R)dV′dV
其中， ∫ₜ∫ₛ∇∇G(R)dV′dV = ∫ₜ∑ᵢn̂ᵢ(∫ᵢR̂(jk + 1/R)G(R)dS′)dV
计算得到结果为并矢::
jη₀/k ∫∫ (k²I + ∇∇)G(R) dV'dV
Kᵣⁿ  =   ∫ Rⁿ dV'
K̂ᵣⁿ  =   ∫ R̂Rⁿ dV'
注意为与两两作用不同，此处加上了 κ 项，因此后续填充时不需加上。
"""
function EFIEOnTetraPWC(tetrat::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    # 场源求积点
    rgt     =   getGQPTetra(tetrat)
    # 采样点对应的体元
    dVt  =   tetrat.volume
    # 常数项
    Jη_0divKdVt =   Params.Jη_0divK * dVt
    ## 计算矩阵元并矢
    # 对源求积点循环
    @inbounds for gi in 1:GQPNTetra
        # 源高斯求积点
        rgi  =      view(rgt, :, gi)
        # 计算 L 算子并矢并乘以权重
        Zts .+=     TetraGQInfo.weight[gi] .* volumeSingularityLOpDyad(rgi, tetrat)
    end

    # 补上常数项
    Zts .*= Jη_0divKdVt*tetrat.κ

    # 计算矩阵对角线项，只在场源重合时出现
    selfImp  =  Params.divJω/tetrat.ε*tetrat.volume
    for ni in 1:3
        Zts[ni, ni] += selfImp
    end

    return Zts
end

"""
采用 PWC 基函数
计算四面体上相关的9个阻抗矩阵元，
此函数方法用于计算场源四面体重合的情况，因此输入有一个四面体信息类型实例
输入：
tetrat   TetrahedraInfo, 场四面体和四面体
计算：
jkη₀∫ₜ∫ₛ(I + 1/k²∇∇)G(R)dV′dV
其中， ∫ₜ∫ₛ∇∇G(R)dV′dV = ∫ₜ∑ᵢn̂ᵢ(∫ᵢR̂(jk + 1/R)G(R)dS′)dV
计算得到结果为并矢::
jη₀/k ∫∫ (k²I + ∇∇)G(R) dV'dV
Kᵣⁿ  =   ∫ Rⁿ dV'
K̂ᵣⁿ  =   ∫ R̂Rⁿ dV'
注意为与两两作用不同，函数将主值积分分开返回以便它用
"""
function EFIEOnTetraPWCSepPV(tetrat::TetrahedraInfo{IT, FT, CT}) where {IT<: Integer, FT<:AbstractFloat, CT<:Complex{FT}}
    # 保存结果的 (3*3) 小数组
    Zts     =   zeros(CT, 3, 3)
    # 场源求积点
    rgt     =   getGQPTetra(tetrat)
    # 采样点对应的体元
    dVt  =   tetrat.volume
    # 常数项
    Jη_0divKdVt =   Params.Jη_0divK * dVt
    ## 计算矩阵元并矢
    # 对源求积点循环
    @inbounds for gi in 1:GQPNTetra
        # 源高斯求积点
        rgi  =      view(rgt, :, gi)
        # 计算 L 算子并矢并乘以权重
        Zts .+=     TetraGQInfo.weight[gi] .* volumeSingularityLOpDyad(rgi, tetrat)
    end

    # 补上常数项
    Zts .*= Jη_0divKdVt

    # 保存主值积分的值，矩阵对角线项，只在场源重合时出现
    ZtsPV   =   Params.divJω*tetrat.volume

    return Zts, ZtsPV
end


"""
本函数用于计算介质体的 PWC 基函数下的 EFIE 阻抗矩阵。
输入信息：
tetrasInfo  :  为包含四面体信息实例的向量
nrwg        :  基函数数目
返回值
Zmat        :  阻抗矩阵
"""
function impedancemat4VIE(tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, nPWC::Integer, bfT::Type{BFT}) where {IT, FT, CT, BFT<:PWC}
    
    # 初始化阻抗矩阵
    Zmat    =   zeros(CT, (nPWC, nPWC))
    # 计算
    impedancemat4VIE!(Zmat, tetrasInfo, bfT)

    return Zmat
    
end


"""
本函数用于计算介质体的 PWC 基函数下的 EFIE 阻抗矩阵。
输入信息：
tetrasInfo  :  为包含四面体信息实例的向量
nrwg        :  基函数数目
返回值
Zmat        :  阻抗矩阵
"""
function impedancemat4VIE!(Zmat::Matrix{CT}, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, ::Type{BFT}; 
    discreteVar = SimulationParams.discreteVar) where {IT, FT, CT, BFT<:PWC}
    
    # 四面体数
    tetrasnum    =   length(tetrasInfo)
    # 几何信息索引区间
    geoInterval =   getGeosInterval(tetrasInfo)
    # 判断体电流的离散方式，
    discreteJ::Bool = discreteVar == "J"

    # 常数
    Rsglr       =   Params.Rsglr
    # Progress Meter
    nbf     =   size(Zmat, 1)
    pmeter  =   Progress(tetrasnum; desc = "Calculating Z (PWC)($nbf × $nbf)")
    # 外层定义为场基函数循环
    @threads for ti in geoInterval
        # 局域的场四面体
        @inbounds local tetrat  =   tetrasInfo[ti]
        # 场四面体介质对比度
        κₜ  =   tetrat.κ

        Rsglrlc =   Rsglr/sqrt(norm(tetrat.ε)/ε_0)
        @inbounds for sj in ti:geoInterval.stop
            # 局域的源三角形
            local tetras  =   tetrasInfo[sj]
            # 源四面体介质对比度
            κₛ  =   tetras.κ
            # 场源距离
            local Rts   =   dist(tetrat.center, tetras.center)

            # isapprox(Rts, Rsglrlc, rtol = 1e-2) && @show ti, sj
            # 判断二者远近，调用不同精度的矩阵元处理函数
            if ti == sj
                Zts, ZtsPV  =   EFIEOnTetraPWCSepPV(tetrat)
                for ni in 1:3
                    n = tetras.inBfsID[ni]
                    for mi in 1:3
                        # 基函数id
                        m = tetrat.inBfsID[mi]
                        # 写入
                        if discreteJ
                            Zmat[m, n]  =   Zts[mi, ni]
                        else
                            Zmat[m, n]  =   Zts[mi, ni]*κₜ
                        end
                    end
                    if discreteJ
                        Zmat[n, n] += ZtsPV/(tetrat.ε - ε_0)
                    else
                        Zmat[n, n] += ZtsPV/tetrat.ε
                    end
                end
            elseif Rts < Rsglrlc
                # 需要进行近奇异性处理的场源四面体
                Zts    =   EFIEOnNearTetrasPWC(tetrat, tetras)
                # 写入数据，利用对称性快速填充，因此要避免重合时重复填充
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m   =   tetrat.inBfsID[mi]
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
                # 计算四面体相关的(3*3)个矩阵元的结果
                Zts    =   EFIEOnTetrasPWC(tetrat, tetras)
                
                # 写入数据
                for ni in 1:3, mi in 1:3
                    # 基函数id
                    m   =   tetrat.inBfsID[mi]
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


# """
# 本函数用于计算介质体的 PWC 基函数下的 EFIE 阻抗矩阵。
# 输入信息：
# tetrasInfo  :  为包含四面体信息实例的向量
# nrwg        :  基函数数目
# 返回值
# Zmat        :  阻抗矩阵
# """
# function impedancemat4VIE!(Zmat::Matrix{CT}, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, ::Type{BFT}) where {IT, FT, CT, BFT<:PWC}
    
#     # 四面体数
#     tetrasnum    =   length(tetrasInfo)
#     isoffset    =   isa(tetrasInfo, OffsetVector)
#     geoInterval = begin 
#         isoffset ? begin
#             st  =   (eachindex(tetrasInfo).offset) + 1
#             st:(st - 1 + tetrasnum)
#         end : begin
#             1:tetrasnum
#         end
#     end
#     # 常数
#     Rsglr       =   Params.Rsglr
#     # Progress Meter
#     nbf     =   size(Zmat, 1)
#     pmeter  =   Progress(tetrasnum; desc = "Calculating Z (PWC)($nbf × $nbf)")
#     # 外层定义为场基函数循环
#     @threads for ti in geoInterval
#         # 局域的场四面体
#         @inbounds local tetrat  =   tetrasInfo[ti]
#         # 场四面体介质对比度
#         κₜ  =   tetrat.κ

#         Rsglrlc =   Rsglr/sqrt(norm(tetrat.ε)/ε_0)
#         @inbounds for sj in ti:geoInterval.stop
#             # 局域的源三角形
#             local tetras  =   tetrasInfo[sj]
#             # 源四面体介质对比度
#             κₛ  =   tetras.κ
#             # 场源距离
#             local Rts   =   dist(tetrat.center, tetras.center)

#             # isapprox(Rts, Rsglrlc, rtol = 1e-2) && @show ti, sj
#             # 判断二者远近，调用不同精度的矩阵元处理函数
#             if ti == sj
#                 Zts    =   EFIEOnTetraPWC(tetrat)
#                 for ni in 1:3, mi in 1:3
#                     # 基函数id
#                     m = tetrat.inBfsID[mi]
#                     n = tetras.inBfsID[ni]
#                     # 写入
#                     Zmat[m, n]  =   Zts[mi, ni]
#                 end
#             elseif Rts < Rsglrlc
#                 # 需要进行近奇异性处理的场源四面体
#                 Zts    =   EFIEOnNearTetrasPWC(tetrat, tetras)
#                 # 写入数据，利用对称性快速填充，因此要避免重合时重复填充
#                 for ni in 1:3, mi in 1:3
#                     # 基函数id
#                     m   =   tetrat.inBfsID[mi]
#                     n   =   tetras.inBfsID[ni]
#                     # 写入
#                     Zmat[m, n]  =   Zts[mi, ni]*κₛ
#                     Zmat[n, m]  =   Zts[mi, ni]*κₜ
#                 end
#             else
#                 # 正常高斯求积
#                 # 计算四面体相关的(3*3)个矩阵元的结果
#                 Zts    =   EFIEOnTetrasPWC(tetrat, tetras)
                
#                 # 写入数据
#                 for ni in 1:3, mi in 1:3
#                     # 基函数id
#                     m   =   tetrat.inBfsID[mi]
#                     n   =   tetras.inBfsID[ni]
#                     # 写入
#                     Zmat[m, n]  =   Zts[mi, ni]*κₛ
#                     Zmat[n, m]  =   Zts[mi, ni]*κₜ
#                 end

#             end # if

#         end #for sj

#         next!(pmeter)

#     end #for ti

#     return nothing
    
# end