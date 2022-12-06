"""
计算三角形重合时的奇异性F1项，即
        ∫∫(1/R)dSdS
的解析值。
输入值：
a, b, c 三角形的三边长
"""
function singularF1(a::FT, b::FT, c::FT) where{FT<:AbstractFloat}
    s   = (a + b + c)/2
    return  -4* (   1/a*log(1 - a/s) + 
                    1/b*log(1 - b/s) + 
                    1/c*log(1 - c/s)    )/3
end

"""
计算三角形重合时的奇异性F1项，即
        ∫∫(1/R)dSdS
的解析值。
输入值：
a, b, c, d 四边形的四边长
"""
function singularF1(a::FT, b::FT, c::FT, d::FT) where{FT<:AbstractFloat}
    s   = (a + b + c + d)/2
    return  -4* (   1/a*log(1 - a/s) + 
                    1/b*log(1 - b/s) + 
                    1/c*log(1 - c/s) +
                    1/c*log(1 - d/s)  )/3
end

"""
计算三角形重合时的奇异性F2项，即
        ∫∫(ρ_m⋅ρ_n/R)dSdS
的解析值。
该函数处理 m==n 即基函数自作用的情况。
输入值：
a, b, c 三角形的三边长
area2   面积的平方
"""
function singularF21(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}
    a2 = a^2; b2 = b^2; c2 = c^2
    s::FT   = (a + b + c)/2
    return  (  (10 - 3*(a2-b2)/c2 - 3*(a2-c2)/b2)*a - 
                (5 - 3*(a2-b2)/c2 - 2*(b2-c2)/a2)*b -
                (5 - 3*(a2-c2)/b2 - 2*(c2-b2)/a2)*c +
                (a2- 3*b2 - 3*c2 - 8*area2/a2)*2/a*log(1 - a/s) +
                (a2- 2*b2 - 4*c2 + 6*area2/b2)*4/b*log(1 - b/s) + 
                (a2- 4*b2 - 2*c2 + 6*area2/c2)*4/c*log(1 - c/s)  )/30
end

"""
计算三角形重合时的奇异性F2项，即
        ∫∫(ρ_m⋅ρ_n/R)dSdS
的解析值。
该函数处理 m!=n 即同一三角形的不同基函数作用的情况。
输入值：
a, b, c 三角形的三边长
area2   面积的平方
"""
function singularF22(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}
    a2 = a^2; b2 = b^2; c2 = c^2
    s::FT   = (a + b + c)/2
    return ((-10 - (a2-b2)/c2 -   (a2-c2)/b2)*a +
            (  5 + (a2-b2)/c2 - 6*(b2-c2)/a2)*b +
            (  5 + (a2-c2)/b2 - 6*(c2-b2)/a2)*c + 
            (  2*a2 -   b2 -   c2 + 4*area2/a2)*12/a*log(1 - a/s) +
            (  9*a2 - 3*b2 -   c2 + 4*area2/b2)*2 /b*log(1 - b/s) +
            (  9*a2 -   b2 - 3*c2 + 4*area2/c2)*2 /c*log(1 - c/s)  )/60
end


"""
面上的近奇异性
rgt, 为场三角形的求积点
tris::TriangleInfo{IT, FT}， 源三角形信息
计算得到结果::
IgS  =   ∫ g(R) dS'   =   ∑ₙ₌₀(coeffgreen(n)*IR[n-1])
"""
function faceSingularityIg(rgt::AbstractVector{FT}, tris::TriangleInfo{IT, FT}, area::FT, facen̂::AbstractVector{FT}) where {IT<:Integer, FT<:Real}
    # 复数数据类型
    CT  =   Complex{FT}
    # 预分配 Ilᵣ  ISᵣ数组
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)

    # 结果
    ISg     =   zero(CT)

    # 开始计算
    # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
    # dts     =   facen̂ ⋅ (rgt .- view(tris.vertices, :, 1))
    dts  =   zero(FT)
    for ii in 1:3
        dts += facen̂[ii] * (rgt[ii] - tris.vertices[ii, 1])
    end
    # 第 gi 个投影点
    r0gi    =   rgt .- dts .* facen̂
    
    # 距离的绝对值、平方
    dtsAbs  =   abs(dts)
    dts²    =   dtsAbs^2

    # 预分配临时变量以加速
    p02jvec =   zero(MVec3D{FT})

    # 对三个边循环
    for edgej in 1:3
        # 构成该边的第一个点
        edgeNodei⁻  =    tris.vertices[:, MoM_Basics.EDGEVmINTriVsID[edgej]]
        # 该边边长
        lⱼ          =   abs(tris.edgel[edgej])
        # 每个边的局部坐标下的 lⱼ⁺ \lⱼ⁻
        # lⱼ⁻     =   (edgeNodei⁻ .- r0gi) ⋅ tris.edgev̂[:, edgej]
        lⱼ⁻  =   zero(FT)
        for ii in 1:3
            lⱼ⁻ += (edgeNodei⁻[ii] - r0gi[ii]) * tris.edgev̂[ii, edgej]
        end
        lⱼ⁺     =   lⱼ⁻ + lⱼ
        # 投影点 r0gi 到各个边的垂足的向量 p02jvec
        # p02jvec =   edgeNodei⁻ .- lⱼ⁻ * tris.edgev̂[:, edgej] .- r0gi
        for ii in 1:3
            p02jvec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * tris.edgev̂[ii, edgej] - r0gi[ii]
        end
        # 该向量长的平方、该向量长
        p02jl   =   p02jvec ⋅ tris.edgen̂[:, edgej]
        p02jl²  =   p02jl^2
        # 将p02ivec化为单位向量
        p02jvec ./=   p02jl
        # 场点到在该边上投影点的距离
        R0²     =   p02jl² + dts²
        # 场点到在该边上正、负端点的距离
        R⁺      =   sqrt(lⱼ⁺^2 + R0²)
        R⁻      =   sqrt(lⱼ⁻^2 + R0²)
        # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
        # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
        let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = 1e-2lⱼ
            if abs(p02jl) < ϵl
                if dtsAbs < ϵl
                    # 视为积分点与边重合，此时 βⱼ = 0
                    continue
                else
                    # 视为投影点与边重合，此时 βⱼ = 0
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                end
            else
                fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                if dtsAbs < ϵl
                    # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                    βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ
                else
                    # 正常处理
                    βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                end
            end #if p02jl

            # 计算 Ilᵣ 所有项
            Ilᵣ[-1]  =   fⱼ
            Ilᵣ[0]   =   lⱼ
            R⁺ⁿ =   one(FT)
            R⁻ⁿ =   one(FT) 
            for n in 1:(SglrOrder-2)
                R⁺ⁿ    *=   R⁺
                R⁻ⁿ    *=   R⁻
                Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
            end # n
            # 计算 ISᵣ 中与边相关的项
            for n in 1:(SglrOrder-2)
                ISᵣ[n]  +=   p02jl * Ilᵣ[n]
            end # n
        end #let
    end # edgej
    # 计算 ISᵣ[n] 与 ISᵣ[n-2] 相关项
    ISᵣ[0]   =   area
    for n in 1:(SglrOrder-2)
        ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
        ISᵣ[n]  /=  n + 2
    end
    # 计算 ISg
    for n in 0:(SglrOrder-1)
        ISg +=   SSCg[n]*ISᵣ[n-1]
    end

    return ISg

end

"""
面上的近奇异性
rgt, 为场三角形的求积点
tris::TriangleInfo{IT, FT}， 源三角形信息
计算得到结果::
IgS  =   ∫ g(R) dS'   =   ∑ₙ₌₀(coeffgreen(n)*IR[n-1])
"""
function faceSingularityIg(rgt::AbstractVector{FT}, tris::Tris4Tetra{IT, FT}, area::FT, facen̂::AbstractVector{FT}) where {IT<:Integer, FT<:Real}
    # 复数数据类型
    CT  =   Complex{FT}
    # 预分配 Ilᵣ  ISᵣ数组
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)

    # 结果
    ISg     =   zero(CT)

    # 开始计算
    # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
    # dts     =   facen̂ ⋅ (rgt .- view(tris.vertices, :, 1))
    dts  =   zero(FT)
    for ii in 1:3
        dts += facen̂[ii] * (rgt[ii] - tris.vertices[ii, 1])
    end
    # 第 gi 个投影点
    r0gi    =   rgt .- dts .* facen̂
    
    # 距离的绝对值、平方
    dtsAbs  =   abs(dts)
    dts²    =   dtsAbs^2
    # 预分配临时变量以加速
    p02jvec =   zero(MVec3D{FT})

    # 对三个边循环
    for edgej in 1:3
        # 构成该边的第一个点
        edgeNodei⁻  =    tris.vertices[:, MoM_Basics.EDGEVmINTriVsID[edgej]]
        # 该边边长
        lⱼ          =   abs(tris.edgel[edgej])
        # 每个边的局部坐标下的 lⱼ⁺ \lⱼ⁻
        # lⱼ⁻     =   (edgeNodei⁻ .- r0gi) ⋅ tris.edgev̂[:, edgej]
        lⱼ⁻  =   zero(FT)
        for ii in 1:3
            lⱼ⁻ += (edgeNodei⁻[ii] - r0gi[ii]) * tris.edgev̂[ii, edgej]
        end
        lⱼ⁺     =   lⱼ⁻ + lⱼ
        # 投影点 r0gi 到各个边的垂足的向量 p02jvec
        # p02jvec =   edgeNodei⁻ .- lⱼ⁻ * tris.edgev̂[:, edgej] .- r0gi
        for ii in 1:3
            p02jvec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * tris.edgev̂[ii, edgej] - r0gi[ii]
        end
        # 该向量长的平方、该向量长
        p02jl   =   p02jvec ⋅ tris.edgen̂[:, edgej]
        p02jl²  =   p02jl^2
        # 将p02ivec化为单位向量
        p02jvec ./=   p02jl
        # 场点到在该边上投影点的距离
        R0²     =   p02jl² + dts²
        # 场点到在该边上正、负端点的距离
        R⁺      =   sqrt(lⱼ⁺^2 + R0²)
        R⁻      =   sqrt(lⱼ⁻^2 + R0²)
        # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
        # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
        let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = 1e-2lⱼ
            if abs(p02jl) < ϵl
                if dtsAbs < ϵl
                    # 视为积分点与边重合，此时 βⱼ = 0
                    continue
                else
                    # 视为投影点与边重合，此时 βⱼ = 0
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                end
            else
                fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                if dtsAbs < ϵl
                    # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                    βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ
                else
                    # 正常处理
                    βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                end
            end #if p02jl

            # 计算 Ilᵣ 所有项
            Ilᵣ[-1]  =   fⱼ
            Ilᵣ[0]   =   lⱼ
            R⁺ⁿ =   one(FT)
            R⁻ⁿ =   one(FT) 
            for n in 1:(SglrOrder-2)
                R⁺ⁿ    *=   R⁺
                R⁻ⁿ    *=   R⁻
                Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
            end # n
            # 计算 ISᵣ 中与边相关的项
            for n in 1:(SglrOrder-2)
                ISᵣ[n]  +=   p02jl * Ilᵣ[n]
            end # n
        end #let
    end # edgej
    # 计算 ISᵣ[n] 与 ISᵣ[n-2] 相关项
    ISᵣ[0]   =   area
    for n in 1:(SglrOrder-2)
        ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
        ISᵣ[n]  /=  n + 2
    end
    # 计算 ISg
    for n in 0:(SglrOrder-1)
        ISg +=   SSCg[n]*ISᵣ[n-1]
    end

    return ISg

end


"""
面上的近奇异性
rgt, 为场求积点
tris::TriangleInfo{IT, FT}， 源多边形信息
计算得到结果::
IgS  =   ∫ g(R) dS'   =   ∑ₙ₌₀(coeffgreen(n)*IR[n-1])
"""
function faceSingularityIg(rgt::AbstractVector{FT}, tris::ST, area::FT, 
    facen̂::AbstractVector{FT}) where {IT<:Integer, FT<:Real, ST<:SurfaceCellType{IT, FT}}
    # 复数数据类型
    CT  =   Complex{FT}
    # 预分配 Ilᵣ  ISᵣ数组
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)

    # 结果
    ISg     =   zero(CT)

    # 开始计算
    # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
    # dts     =   facen̂ ⋅ (rgt .- view(tris.vertices, :, 1))
    dts  =   zero(FT)
    for ii in 1:3
        dts += facen̂[ii] * (rgt[ii] - tris.vertices[ii, 1])
    end
    # 第 gi 个投影点
    r0gi    =   rgt .- dts .* facen̂
    
    # 距离的绝对值、平方
    dtsAbs  =   abs(dts)
    dts²    =   dtsAbs^2
    # 边数
    edgeN   =   length(tris.edgel)
    
    # 预分配临时变量以加速
    p02jvec =   zero(MVec3D{FT})

    # 对三个边循环
    for edgej in 1:edgeN
        # 构成该边的第一个点
        edgeNodei⁻  =    tris.vertices[:, edgej]
        # 该边边长
        lⱼ          =   abs(tris.edgel[edgej])
        # 每个边的局部坐标下的 lⱼ⁺ \lⱼ⁻
        # lⱼ⁻     =   (edgeNodei⁻ .- r0gi) ⋅ tris.edgev̂[:, edgej]
        lⱼ⁻  =   zero(FT)
        for ii in 1:3
            lⱼ⁻ += (edgeNodei⁻[ii] - r0gi[ii]) * tris.edgev̂[ii, edgej]
        end
        lⱼ⁺     =   lⱼ⁻ + lⱼ
        # 投影点 r0gi 到各个边的垂足的向量 p02jvec
        # p02jvec =   edgeNodei⁻ .- lⱼ⁻ * tris.edgev̂[:, edgej] .- r0gi
        for ii in 1:3
            p02jvec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * tris.edgev̂[ii, edgej] - r0gi[ii]
        end
        # 该向量长的平方、该向量长
        p02jl   =   p02jvec ⋅ tris.edgen̂[:, edgej]
        p02jl²  =   p02jl^2
        # 将p02ivec化为单位向量
        p02jvec ./=   p02jl
        # 场点到在该边上投影点的距离
        R0²     =   p02jl² + dts²
        # 场点到在该边上正、负端点的距离
        R⁺      =   sqrt(lⱼ⁺^2 + R0²)
        R⁻      =   sqrt(lⱼ⁻^2 + R0²)
        # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
        # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
        let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = 1e-2lⱼ
            if abs(p02jl) < ϵl
                if dtsAbs < ϵl
                    # 视为积分点与边重合，此时 βⱼ = 0
                    continue
                else
                    # 视为投影点与边重合，此时 βⱼ = 0
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                end
            else
                fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                if dtsAbs < ϵl
                    # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                    βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ
                else
                    # 正常处理
                    βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                end
            end #if p02jl

            # 计算 Ilᵣ 所有项
            Ilᵣ[-1]  =   fⱼ
            Ilᵣ[0]   =   lⱼ
            R⁺ⁿ =   one(FT)
            R⁻ⁿ =   one(FT) 
            for n in 1:(SglrOrder-2)
                R⁺ⁿ    *=   R⁺
                R⁻ⁿ    *=   R⁻
                Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
            end # n
            # 计算 ISᵣ 中与边相关的项
            for n in 1:(SglrOrder-2)
                ISᵣ[n]  +=   p02jl * Ilᵣ[n]
            end # n
        end #let
    end # edgej
    # 计算 ISᵣ[n] 与 ISᵣ[n-2] 相关项
    ISᵣ[0]   =   area
    for n in 1:(SglrOrder-2)
        ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
        ISᵣ[n]  /=  n + 2
    end
    # 计算 ISg
    for n in 0:(SglrOrder-1)
        ISg +=   SSCg[n]*ISᵣ[n-1]
    end

    return ISg

end


"""
面上的近奇异性
rgt, 为场三角形的求积点
tris::TriangleInfo{IT, FT}， 源三角形信息
计算得到结果::
IgS     =   ∫ g(R) dS'      =   ∑ₙ₌₀(coeffgreen(n)*IR[n-1])
IvecgS  =   ∫ Rvec g(R) dS' =   ∑ₗⱼ ûⱼ ∑ₙ₌₀(coeffgreen(n)/(n+1)*Ilᵣ[n-1]) + dn̂IgS
"""
function faceSingularityIgIvecg(rgt::AbstractVector{FT}, tris::ST, area::FT, 
    facen̂::AbstractVector{FT}) where {IT<:Integer, FT<:Real, ST<:SurfaceCellType{IT, FT}}

    # 复数数据类型
    CT  =   Complex{FT}
    # 预分配 Ilᵣ  ISᵣ数组
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)

    # 结果
    ISg     =   zero(CT)
    IvecSg  =   zero(MVec3D{CT})

    # 开始计算
    # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
    # dts     =   facen̂ ⋅ (rgt .- view(tris.vertices, :, 1))
    dts  =   zero(FT)
    for ii in 1:3
        dts += facen̂[ii] * (rgt[ii] - tris.vertices[ii, 1])
    end
    # 第 gi 个投影点
    r0gi    =   rgt .- dts .* facen̂
    
    # 距离的绝对值、平方
    dtsAbs  =   abs(dts)
    dts²    =   dtsAbs^2
    # 预分配临时变量以加速
    p02jvec =   zero(MVec3D{FT})

    # 对三个边循环
    for edgej in 1:3
        # 构成该边的第一个点
        edgeNodei⁻  =    tris.vertices[:, MoM_Basics.EDGEVmINTriVsID[edgej]]
        # 该边边长
        lⱼ          =   abs(tris.edgel[edgej])
        # 每个边的局部坐标下的 lⱼ⁺ \lⱼ⁻
        # lⱼ⁻     =   (edgeNodei⁻ .- r0gi) ⋅ tris.edgev̂[:, edgej]
        lⱼ⁻  =   zero(FT)
        for ii in 1:3
            lⱼ⁻ += (edgeNodei⁻[ii] - r0gi[ii]) * tris.edgev̂[ii, edgej]
        end
        lⱼ⁺     =   lⱼ⁻ + lⱼ
        # 投影点 r0gi 到各个边的垂足的向量 p02jvec
        # p02jvec =   edgeNodei⁻ .- lⱼ⁻ * tris.edgev̂[:, edgej] .- r0gi
        for ii in 1:3
            p02jvec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * tris.edgev̂[ii, edgej] - r0gi[ii]
        end
        # 该向量长的平方、该向量长
        p02jl   =   p02jvec ⋅ tris.edgen̂[:, edgej]
        p02jl²  =   p02jl^2
        # 将p02ivec化为单位向量
        p02jvec ./=   p02jl
        # 场点到在该边上投影点的距离
        R0²     =   p02jl² + dts²
        # 场点到在该边上正、负端点的距离
        R⁺      =   sqrt(lⱼ⁺^2 + R0²)
        R⁻      =   sqrt(lⱼ⁻^2 + R0²)
        # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
        # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
        let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = 1e-2lⱼ
            if abs(p02jl) < ϵl
                if dtsAbs < ϵl
                    # 视为积分点与边重合，此时 βⱼ = 0
                    continue
                else
                    # 视为投影点与边重合，此时 βⱼ = 0
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                end
            else
                fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                if dtsAbs < ϵl
                    # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                    βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ
                else
                    # 正常处理
                    βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                end
            end #if p02jl

            # 计算 Ilᵣ 所有项
            Ilᵣ[-1]  =   fⱼ
            Ilᵣ[0]   =   lⱼ
            R⁺ⁿ =   one(FT)
            R⁻ⁿ =   one(FT) 
            for n in 1:(SglrOrder-2)
                R⁺ⁿ    *=   R⁺
                R⁻ⁿ    *=   R⁻
                Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
            end # n
            # 计算 ISᵣ 中与边相关的项
            for n in 1:(SglrOrder-2)
                ISᵣ[n]  +=   p02jl * Ilᵣ[n]
            end # n
            # 计算 ∑ₙ₌₀(coeffgreen(n)/(n+1)*Ilᵣ[n-1])
            Cdvnp1I     =    CT(-Ilᵣ[1])
            for n in 1:(SglrOrder-3)
                Cdvnp1I -=  SSCgdivnp1[n] * Ilᵣ[n + 1]
            end

            # 累加 IvecSg 与边有关的项
            @views IvecSg   .+=  tris.edgen̂[:, edgej] .* Cdvnp1I
        end #let


    end # edgej
    # 计算 ISᵣ[n] 与 ISᵣ[n-2] 相关项
    ISᵣ[0]   =   area
    for n in 1:(SglrOrder-2)
        ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
        ISᵣ[n]  /=  n + 2
    end
    # 计算 ISg
    for n in 0:(SglrOrder-1)
        ISg +=   SSCg[n]*ISᵣ[n-1]
    end

    # 累加 IvecSg 与面有关的项
    IvecSg   .+=  facen̂ * (dts*ISg)

    return ISg, IvecSg

end

"""
面上的近奇异性
rgt, 为场三角形的求积点
tris::TriangleInfo{IT, FT}， 源三角形信息
计算得到结果::
IgS     =   ∫ g(R) dS'      =   ∑ₙ₌₀(coeffgreen(n)*IR[n-1])
IvecgS  =   ∫ Rvec g(R) dS' =   ∑ₗⱼ ûⱼ ∑ₙ₌₀(coeffgreen(n)/(n+1)*Ilᵣ[n-1]) + dn̂IgS
I∇gS    =   ∫ ∇g(R) dS'     =   ∑ₙ₌₀VSC₃ⁿ*(-1/(n+2)∑ₗⱼûⱼIlᵣ[n+2] + dn̂IgS )
"""
function faceSingularityIgIvecgI∇gS(rgt::AbstractVector{FT}, tris::ST, area::FT, 
    facen̂::AbstractVector{FT}) where {IT<:Integer, FT<:Real, ST<:SurfaceCellType{IT, FT}}

    # 复数数据类型
    CT  =   Complex{FT}
    # 预分配 Ilᵣ  ISᵣ数组
    Ilᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)
    ISᵣ     =   OffsetArray(zero(MVector{SglrOrder, FT}), -2)

    # 结果
    ISg     =   zero(CT)
    IvecSg  =   zero(MVec3D{CT})
    I∇gS    =   zero(MVec3D{CT})

    # 开始计算
    # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
    # dts     =   facen̂ ⋅ (rgt .- view(tris.vertices, :, 1))
    dts  =   zero(FT)
    for ii in 1:3
        dts += facen̂[ii] * (rgt[ii] - tris.vertices[ii, 1])
    end
    # 第 gi 个投影点
    r0gi    =   rgt .- dts .* facen̂
    
    # 距离的绝对值、平方
    dtsAbs  =   abs(dts)
    dts²    =   dtsAbs^2

    # 预分配临时变量以加速
    p02jvec =   zero(MVec3D{FT})

    # 对三个边循环
    for edgej in 1:3
        # 构成该边的第一个点
        edgeNodei⁻  =    tris.vertices[:, MoM_Basics.EDGEVmINTriVsID[edgej]]
        # 该边边长
        lⱼ          =   abs(tris.edgel[edgej])
        # 每个边的局部坐标下的 lⱼ⁺ \lⱼ⁻
        # lⱼ⁻     =   (edgeNodei⁻ .- r0gi) ⋅ tris.edgev̂[:, edgej]
        lⱼ⁻  =   zero(FT)
        for ii in 1:3
            lⱼ⁻ += (edgeNodei⁻[ii] - r0gi[ii]) * tris.edgev̂[ii, edgej]
        end
        lⱼ⁺     =   lⱼ⁻ + lⱼ
        # 投影点 r0gi 到各个边的垂足的向量 p02jvec
        # p02jvec =   edgeNodei⁻ .- lⱼ⁻ * tris.edgev̂[:, edgej] .- r0gi
        for ii in 1:3
            p02jvec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * tris.edgev̂[ii, edgej] - r0gi[ii]
        end
        # 该向量长的平方、该向量长
        p02jl   =   p02jvec ⋅ tris.edgen̂[:, edgej]
        p02jl²  =   p02jl^2
        # 将p02ivec化为单位向量
        p02jvec ./=   p02jl
        # 场点到在该边上投影点的距离
        R0²     =   p02jl² + dts²
        # 场点到在该边上正、负端点的距离
        R⁺      =   sqrt(lⱼ⁺^2 + R0²)
        R⁻      =   sqrt(lⱼ⁻^2 + R0²)
        # 投影点过于靠近该边即 p02jl = 0 时， βⱼ = 0
        # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
        let fⱼ = zero(FT), βⱼ =  zero(FT), ϵl = 1e-2lⱼ
            if abs(p02jl) < ϵl
                if dtsAbs < ϵl
                    # 视为积分点与边重合，此时 βⱼ = 0
                    continue
                else
                    # 视为投影点与边重合，此时 βⱼ = 0
                    fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                end
            else
                fⱼ      =   log((lⱼ⁺ + R⁺)/(lⱼ⁻ + R⁻))
                if dtsAbs < ϵl
                    # 视为积分点与源三角形面重合，此时 βⱼ ≠ 0
                    βⱼ      =   atan((p02jl*lⱼ⁺)/R0²) - atan((p02jl*lⱼ⁻)/R0²)
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ
                else
                    # 正常处理
                    βⱼ      =   atan((p02jl*lⱼ⁺)/(R0² + dtsAbs*R⁺)) - atan((p02jl*lⱼ⁻)/(R0² + dtsAbs*R⁻))
                    # 累加 ISᵣ[-1] 项
                    ISᵣ[-1] +=   p02jl*fⱼ - dtsAbs*βⱼ
                end
            end #if p02jl

            # 计算 Ilᵣ 所有项
            Ilᵣ[-1]  =   fⱼ
            Ilᵣ[0]   =   lⱼ
            R⁺ⁿ =   one(FT)
            R⁻ⁿ =   one(FT) 
            for n in 1:(SglrOrder-2)
                R⁺ⁿ    *=   R⁺
                R⁻ⁿ    *=   R⁻
                Ilᵣ[n]  =   (lⱼ⁺*R⁺ⁿ - lⱼ⁻*R⁻ⁿ  + n*R0²*Ilᵣ[n-2])/(n+1)
            end # n
            # 计算 ISᵣ 中与边相关的项
            for n in 1:(SglrOrder-2)
                ISᵣ[n]  +=   p02jl * Ilᵣ[n]
            end # n
            # 计算 ∑ₙ₌₀(coeffgreen(n)/(n+1)*Ilᵣ[n-1])
            Cdvnp1I     =    CT(-Ilᵣ[1])
            for n in 1:(SglrOrder-3)
                Cdvnp1I -=  SSCgdivnp1[n] * Ilᵣ[n + 1]
            end

            # 累加 IvecSg 与边有关的项
            @views IvecSg   .+=  tris.edgen̂[:, edgej] .* Cdvnp1I
        end #let
        C2I =  Ilᵣ ⋅ VSC₃ⁿ

        # 累加 ûC₂Ipn̂dC₂Kvec 与边有关的项
        @views I∇gS   .+=  edgen̂[:, edgei] .* C2I
    end # edgej

    # 其它 ISᵣ 相关项
    diKᵣm3  *=   sign(dts)
    ISᵣ[0]   =   abs(volumeCell.facesArea[iface])
    for n in 1:(SglrOrder-1)
        ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
        ISᵣ[n]  /=  n + 2
    end
    #= 累加面 格林函数在 面上的 奇异性计算结果, 注意此处点乘不能改变顺序
    Julia 程序处理 点乘时，若前面一个向量为复数，则取其共轭，因此，要将 浮点数 向量 ISᵣ 放在前面 =#
    dC₁K    +=  dts*div4π*(ISᵣ ⋅ VSC₁ⁿ)

    # 累加 ûC₂Ipn̂dC₂Kvec 与面有关的项
    @views I∇gS    .+=  (diKᵣm3 + dts*(ISᵣ[-1:(SglrOrder-3)] ⋅ VSC₂ⁿ)) .* volumeCell.facesn̂[:, iface]

    # 计算 ISᵣ[n] 与 ISᵣ[n-2] 相关项
    ISᵣ[0]   =   area
    for n in 1:(SglrOrder-2)
        ISᵣ[n]  +=  n*dts²*ISᵣ[n-2]
        ISᵣ[n]  /=  n + 2
    end
    # 计算 ISg
    for n in 0:(SglrOrder-1)
        ISg +=   SSCg[n]*ISᵣ[n-1]
    end

    # 累加 IvecSg 与面有关的项
    IvecSg   .+=  facen̂ * (dts*ISg)

    return ISg, IvecSg, I∇gS

end


"""
面上的近奇异性
rgts::MMatrix{GQPNTriSglr, 3, Complex{FT}}, 为场三角形的所有高斯求积点
tris::TriangleInfo{IT, FT}， 源三角形信息
计算得到结果::
r0tProj2s::MMatrix{3, GQPNTriSglr, Complex{FT}, 3GQPNTriSglr}, 积分点在源三角形上的投影点
Iᵣ  =   ∫ 1/R dS'   =   ∑₁³(p02il*fᵢ - dtsAbs*βᵢ)
Iᵨ  =   ∫ ρ/R dS'   =   0.5∑₁³{ûᵢ[(R0^2*fᵢ + li⁺*R⁺ - li⁻*R⁻)]}
"""
function faceSingularityIᵣIᵨ(rgt::MMatrix{3, GQPNTriSglr, FT, 3GQPNTriSglr}, 
    tris::TriangleInfo{IT, FT}) where {IT<:Integer, FT<:Real}
    # 预分配场积分点到其在源三角形所有投影点
    r0tProj2s   =   zero(MMatrix{3, GQPNTriSglr, FT, 3GQPNTriSglr})
    # 预分配 Iᵣ、Iᵨ 结果数组
    Iᵣ  =   zero(MVector{GQPNTriSglr, FT})
    Iᵨ  =   zero(MMatrix{3, GQPNTriSglr, FT})
    # 预分配临时变量以加速
    p02ivec =   zero(MVec3D{FT})
    
    # 对高斯求积点循环
    for gi in 1:GQPNTriSglr
        # 第gi个求积点
        rgit     =   rgt[:, gi]
        # 该积分点在到源三角形上的距离（带正负，以源三角形法向为正向）
        # dts      =   tris.facen̂ ⋅ (rgit .- tris.vertices[:, 1])
        dts  =   zero(FT)
        for ii in 1:3
            dts += tris.facen̂[ii] * (rgit[ii] - tris.vertices[ii, 1])
        end
        # 第 gi 个投影点
        r0gi     =   view(r0tProj2s, :, gi)
        r0gi    .=   rgit .- dts .* tris.facen̂
        
        # 距离的绝对值、平方
        dtsAbs  =   abs(dts)
        dts²    =   dtsAbs^2

        # 对三个边循环
        for edgei in 1:3
            # 构成该边的第一个点
            edgeNodei⁻  =    tris.vertices[:, MoM_Basics.EDGEVmINTriVsID[edgei]]
            # 该边边长
            li          =   abs(tris.edgel[edgei])
            # 每个边的局部坐标下的 li⁺ \li⁻
            # li⁻  =   (edgeNodei⁻ .- r0gi) ⋅ tris.edgev̂[:, edgei]
            li⁻  =   zero(FT)
            for ii in 1:3
                li⁻ += (edgeNodei⁻[ii] - r0gi[ii]) * tris.edgev̂[ii, edgei]
            end
            li⁺         =   li⁻ + li
            # 投影点 r0gi 到各个边的垂足的向量 p02ivec
            #  p02ivec =   edgeNodei⁻ .- li⁻ * tris.edgev̂[:, edgei] .- r0gi
            for ii in 1:3
                p02ivec[ii] = edgeNodei⁻[ii] - li⁻ * tris.edgev̂[ii, edgei] - r0gi[ii]
            end
            # 该向量长的平方、该向量长
            p02il   =   p02ivec ⋅ tris.edgen̂[:, edgei]
            p02il²  =   p02il^2
            # 将p02ivec化为单位向量
            p02ivec ./=   p02il
            # 场点到在该边上投影点的距离
            R0²     =   p02il² + dts²
            # 场点到在该边上正、负端点的距离
            R⁺      =   sqrt(li⁺^2 + R0²)
            R⁻      =   sqrt(li⁻^2 + R0²)
            # 投影点过于靠近该边即 p02il = 0 时， βᵢ = 0
            # 此处为避免数值误差，将阈值 ϵl 设定为 1e-2 边长
            let fᵢ = zero(FT), βᵢ =  zero(FT), ϵl = 1e-3li
                if abs(p02il) < ϵl
                    if dtsAbs < ϵl
                        # 视为积分点与边重合，此时 βᵢ = 0， Iᵣ不变，Iᵨ如下
                        Iᵨ[:, gi] .+=    tris.edgen̂[:, edgei] * (0.5(li⁺*R⁺ - li⁻*R⁻))
                    else
                        # 视为投影点与边重合，此时 βᵢ = 0， Iᵣ不变，Iᵨ如下
                        fᵢ      =   log((li⁺ + R⁺)/(li⁻ + R⁻))
                        Iᵨ[:, gi] .+=    tris.edgen̂[:, edgei] * (0.5(R0²*fᵢ + li⁺*R⁺ - li⁻*R⁻))
                    end
                else
                    fᵢ      =   log((li⁺ + R⁺)/(li⁻ + R⁻))
                    if dtsAbs < ϵl
                        # 视为积分点与源三角形面重合，此时 βᵢ ≠ 0
                        βᵢ      =   atan((p02il*li⁺)/R0²) - atan((p02il*li⁻)/R0²)
                        # Iᵣ累加
                        Iᵣ[gi] +=   p02il*fᵢ
                    else
                        # 正常处理
                        βᵢ      =   atan((p02il*li⁺)/(R0² + dtsAbs*R⁺)) - atan((p02il*li⁻)/(R0² + dtsAbs*R⁻))
                        Iᵣ[gi] +=   p02il*fᵢ - dtsAbs*βᵢ
                    end
                    # Iᵨ累加
                    Iᵨ[:, gi] .+=    tris.edgen̂[:, edgei] * (0.5(R0²*fᵢ + li⁺*R⁺ - li⁻*R⁻))
                end #if p02il
            end #let
        end # edgei
    end #gi

    return r0tProj2s, Iᵣ, Iᵨ

end