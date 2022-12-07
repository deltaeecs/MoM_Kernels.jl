"""
在球坐标为r̂θϕ处计算辐射积分，采用RWG基函数时，三角形上没有统一的电流值，每一点上都是三边电流的叠加，
此时:
N(θ, ϕ) =   ∑ₙ(∫ₛ Jˢ exp( jkr̂(θ, ϕ)⋅rₙ ) dS)
        =   ∑ₙ(∫ₛ (∑ₜₙ₌₁³ Iₙfₙ)exp(jkr̂(θ, ϕ)⋅rₙ) dS)
        =   ∑ₙ(Sₜ (∑ₜₙ₌₁³ Iₙlₙρₙ/(2Sₙ))exp(jkr̂(θ, ϕ)⋅rₙ) )
        =   ∑ₙ(∑ᵢWᵢ(∑ₜₙ₌₁³ Iₙlₙρₙ/2)exp(jkr̂(θ, ϕ)⋅rₙ) )
"""
function raditionalIntegralNθϕCal(r̂θϕ::r̂θϕInfo{FT}, trianglesInfo::Vector{TriangleInfo{IT, FT}}, 
                                            Jtris::Array{CT, 3}) where {IT<:Integer, FT<:Real, CT<:Complex}
    # 直角坐标结果
    Nxyz    =   zero(MVec3D{CT})
    # 球坐标θϕ分量结果
    Nθϕ     =   zero(MVector{2, CT})
    # 解压这角度相关信息
    r̂, θhat, ϕhat, θϕ   =   r̂θϕ.r̂, r̂θϕ.θhat, r̂θϕ.ϕhat, r̂θϕ.θϕ
    # 常数
    JK_0 = Params.JK_0
    # 对三角形循环计算
    @inbounds for ti in 1:length(trianglesInfo)
        # 第ti个三角形
        tri     =   trianglesInfo[ti]
        # 该三角形所在的三个基函数id
        # tms     =   tri.inBfsID
        # 相关的电流
        Jtri    =   @view Jtris[:, :, ti]
        # 积分值临时变量
        JSexp   =   zero(MVec3D{CT})

        # 对高斯求积点循环
        for gi in 1:GQPNTri    
            # 采样点
            rgi     =   getGQPTri(tri, gi)
            # 积分值累加
            @views JSexp .+=   Jtri[:, gi] .* (exp(JK_0 * (r̂ ⋅ rgi)) * TriGQInfo.weight[gi])
        end #for gi
        # 结果修正
        JSexp     .*=   tri.area
        # 累加到辐射积分结果上
        Nxyz      .+=   JSexp
    end #for ti

    # 将直角坐标转换到球坐标
    Nθϕ[1]  =   θhat ⋅ Nxyz
    Nθϕ[2]  =   ϕhat ⋅ Nxyz

    return Nθϕ
end

"""
在球坐标为r̂θϕ处计算辐射积分，采用RWG基函数时，三角形上没有统一的电流值，每一点上都是三边电流的叠加，
此时:
N(θ, ϕ) =   ∑ₙ(∫ₛ Jˢ exp( jkr̂(θ, ϕ)⋅rₙ ) dS)
        =   ∑ₙ(∫ₛ (∑ₜₙ₌₁³ Iₙfₙ)exp(jkr̂(θ, ϕ)⋅rₙ) dS)
        =   ∑ₙ(Sₜ (∑ₜₙ₌₁³ Iₙlₙρₙ/(2Sₙ))exp(jkr̂(θ, ϕ)⋅rₙ) )
        =   ∑ₙ(∑ᵢWᵢ(∑ₜₙ₌₁³ Iₙlₙρₙ/2)exp(jkr̂(θ, ϕ)⋅rₙ) )
"""
function raditionalIntegralNθϕCal(r̂θϕ::r̂θϕInfo{FT}, trianglesInfo::Vector{TriangleInfo{IT, FT}}, 
                                            Jtris::Array{CT, 2}) where {IT<:Integer, FT<:Real, CT<:Complex}
    # 直角坐标结果
    Nxyz    =   zero(MVec3D{CT})
    # 球坐标θϕ分量结果
    Nθϕ     =   zero(MVector{2, CT})
    # 解压这角度相关信息
    r̂, θhat, ϕhat, θϕ   =   r̂θϕ.r̂, r̂θϕ.θhat, r̂θϕ.ϕhat, r̂θϕ.θϕ
    # 常数
    JK_0 = Params.JK_0
    # 对三角形循环计算
    @inbounds for ti in 1:length(trianglesInfo)
        # 第ti个三角形
        tri     =   trianglesInfo[ti]
        # 该三角形所在的三个基函数id
        # tms     =   tri.inBfsID
        # 相关的电流
        Jtri    =   @view Jtris[:, ti]
        # 积分值临时变量
        JSexp   =   zero(MVec3D{CT})

        # 对高斯求积点循环
        for gi in 1:GQPNTri    
            # 采样点
            rgi     =   getGQPTri(tri, gi)
            # 积分值累加
            @views JSexp .+=   Jtri .* (exp(JK_0 * (r̂ ⋅ rgi)) * TriGQInfo.weight[gi])
        end #for gi
        # 结果修正
        JSexp     .*=   tri.area
        # 累加到辐射积分结果上
        Nxyz      .+=   JSexp
    end #for ti

    # 将直角坐标转换到球坐标
    Nθϕ[1]  =   θhat ⋅ Nxyz
    Nθϕ[2]  =   ϕhat ⋅ Nxyz

    return Nθϕ
end

"""
在球坐标为r̂θϕ处计算辐射积分，采用 SWG 基函数时，四面体上没有统一的电流值，每一点上都是四面电流的叠加，
此时:
N(θ, ϕ) =   ∑ₙ(∫ₛ Jˢ exp( jkr̂(θ, ϕ)⋅rₙ ) dV)
        =   ∑ₙ(∫ₛ (∑ₜₙ₌₁⁴ Iₙfₙ)exp(jkr̂(θ, ϕ)⋅rₙ) dV)
        =   ∑ₙ(Sₜ (∑ₜₙ₌₁⁴ IₙSₙρₙ/(3Vₙ))exp(jkr̂(θ, ϕ)⋅rₙ) )
        =   ∑ₙ(∑ᵢWᵢ(∑ₜₙ₌₁⁴ IₙSₙρₙ/3)exp(jkr̂(θ, ϕ)⋅rₙ) )
"""
function raditionalIntegralNθϕCal(r̂θϕ::r̂θϕInfo{FT}, tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}, 
                                            Js::Array{CT, 2}) where {IT<:Integer, FT<:Real, CT<:Complex}
    # 直角坐标结果
    Nxyz    =   zero(MVec3D{CT})
    # 球坐标θϕ分量结果
    Nθϕ     =   zero(MVector{2, CT})
    # 解压这角度相关信息
    r̂, θhat, ϕhat, θϕ   =   r̂θϕ.r̂, r̂θϕ.θhat, r̂θϕ.ϕhat, r̂θϕ.θϕ
    # 常数
    JK_0 = Params.JK_0
    # 四面体数
    ntetra   =   length(tetrasInfo)
    # 是否为偏置数组
    isoffset    =   isa(tetrasInfo, OffsetVector)
    geoIdx      =   eachindex(tetrasInfo)
    # 对六面体循环计算
    for it in 1:ntetra
        # 第it个六面体
        idx = isoffset ? (geoIdx.offset + it) : it
        tetra   =   tetrasInfo[idx]
        # 该三角形所在的三个基函数id
        # tms     =   tetra.inBfsID
        # 相关的电流
        Jti     =   @view Js[:, it]
        # 积分值临时变量
        Jtexp   =   zero(MVec3D{CT})

        # 对高斯求积点循环
        for gi in 1:GQPNTetra    
            # 采样点
            rgi     =   getGQPTetra(tetra, gi)
            # 积分值累加
            @views Jtexp .+=   Jti .* (exp(JK_0 * (r̂ ⋅ rgi)) * TetraGQInfo.weight[gi])
        end #for gi
        # 结果修正
        Jtexp   .*=   tetra.volume
        # 累加到辐射积分结果上
        Nxyz    .+=   Jtexp
    end #for it

    # 将直角坐标转换到球坐标
    Nθϕ[1]  =   θhat ⋅ Nxyz
    Nθϕ[2]  =   ϕhat ⋅ Nxyz

    return Nθϕ
end

"""
在球坐标为r̂θϕ处计算辐射积分，采用 RBF 基函数时，六面体上没有统一的电流值
N(θ, ϕ) =   ∑ₙ(∫ₛ Jˢ exp( jkr̂(θ, ϕ)⋅rₙ ) dV)
"""
function raditionalIntegralNθϕCal(r̂θϕ::r̂θϕInfo{FT}, hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}, 
                                            Js::Array{CT, 2}) where {IT<:Integer, FT<:Real, CT<:Complex}
    # 直角坐标结果
    Nxyz    =   zero(MVec3D{CT})
    # 球坐标θϕ分量结果
    Nθϕ     =   zero(MVector{2, CT})
    # 解压这角度相关信息
    r̂, θhat, ϕhat   =   r̂θϕ.r̂, r̂θϕ.θhat, r̂θϕ.ϕhat
    # 常数
    JK_0 = Params.JK_0
    # 四面体数
    nhexa   =   length(hexasInfo)
    # 是否为偏置数组
    isoffset    =   isa(hexasInfo, OffsetVector)
    geoIdx      =   eachindex(hexasInfo)
    # 对六面体循环计算
    for it in 1:nhexa
        # 第it个六面体
        idx = isoffset ? (geoIdx.offset + it) : it
        hexa   =   hexasInfo[idx]
        # 该三角形所在的三个基函数id
        # tms     =   hexa.inBfsID
        # 相关的电流
        Jti     =   @view Js[:, it]
        # 积分值临时变量
        Jtexp   =   zero(MVec3D{CT})

        # 对高斯求积点循环
        for gi in 1:GQPNHexa    
            # 采样点
            rgi     =   getGQPHexa(hexa, gi)
            # 积分值累加
            @views Jtexp .+=   Jti .* (exp(JK_0 * (r̂ ⋅ rgi)) * HexaGQInfo.weight[gi])
        end #for gi
        # 结果修正
        Jtexp   .*=   hexa.volume
        # 累加到辐射积分结果上
        Nxyz    .+=   Jtexp
    end #for it

    # 将直角坐标转换到球坐标
    Nθϕ[1]  =   θhat ⋅ Nxyz
    Nθϕ[2]  =   ϕhat ⋅ Nxyz

    return Nθϕ
end

"""
在设定好的观测角度上的球坐标处计算辐射积分，采用RWG基函数时，三角形上没有统一的电流值，每一点上都是三边电流的叠加，
此时:
N(θ, ϕ) =   ∑ₙ(∫ₛ Jˢ exp( jkr̂(θ, ϕ)⋅rₙ ) dS)
        =   ∑ₙ(∫ₛ (∑ₜₙ₌₁³ Iₙfₙ)exp(jkr̂(θ, ϕ)⋅rₙ) dS)
        =   ∑ₙ(Sₜ (∑ₜₙ₌₁³ Iₙlₙρₙ/(2Sₙ))exp(jkr̂(θ, ϕ)⋅rₙ) )
        =   ∑ₙ(∑ᵢWᵢ(∑ₜₙ₌₁³ Iₙlₙρₙ/2)exp(jkr̂(θ, ϕ)⋅rₙ) )
"""
function raditionalIntegralNCal(θs_obs::LinRange{FT}, ϕs_obs::LinRange{FT},
                                            geosInfo::Vector{TriangleInfo{IT, FT}},
                                            Jgeos::Array{CT}) where {IT<:Integer, FT<:Real, CT<:Complex}


    # 观测角度信息
    Nθ_obs      =   length(θs_obs)
    Nϕ_obs      =   length(ϕs_obs)
    nobs        =   length(θs_obs) * length(ϕs_obs)
    θsobsInfo   =   [∠Info{FT}(θ_obs) for θ_obs in θs_obs]
    ϕsobsInfo   =   [∠Info{FT}(ϕ_obs) for ϕ_obs in ϕs_obs]
    r̂θsϕs       =   [r̂θϕInfo{FT}(θobsInfo, ϕobsInfo) for θobsInfo in θsobsInfo, ϕobsInfo in ϕsobsInfo]
    
    # 预分配辐射积分内存
    Nθsϕs       =   zeros(CT, (2, length(θsobsInfo), length(ϕsobsInfo)))
    Nθsϕsrsp    =   reshape(Nθsϕs, (2, nobs))

    # Progress Meter
    pmeter  =   Progress(nobs, "Calculating radiational integral ($Nθ_obs × $Nϕ_obs)")
    # 输入的为高斯求积点电流则提取中心值
    Jgeos2 = if typeof(Jgeos)<:Array{CT, 3}
        Jgeos[:, 1, :]
    elseif typeof(Jgeos)<:Array{CT, 3}
        Jgeos
    else
        throw("Error in current calculation, check it!")
    end
    # 对观测角度循环计算
    @threads for ii in 1:nobs
        # 调用函数计算
        @inbounds Nθsϕsrsp[:, ii]   .=  raditionalIntegralNθϕCal(r̂θsϕs[ii], geosInfo, Jgeos2)
        next!(pmeter)
    end

    return Nθsϕs

end
