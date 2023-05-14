@doc raw"""
    singularF1(a::FT, b::FT, c::FT) where{FT<:AbstractFloat}

计算边长为`a, b, c`的三角形重合时的奇异性F1项，即
``\int{\int{\frac{1}{R}}dS}``
的解析值。
"""
function singularF1(a::FT, b::FT, c::FT) where{FT<:AbstractFloat}
    s   = (a + b + c)/2
    return  -4* (   1/a*log(1 - a/s) + 
                    1/b*log(1 - b/s) + 
                    1/c*log(1 - c/s)    )/3
end

@doc raw"""
    singularF1(a::FT, b::FT, c::FT, d::FT) where{FT<:AbstractFloat}

计算边长为`a, b, c, d`的四边形重合时的奇异性F1项，即
``\int{\int{\frac{1}{R}}dS}``
的解析值。
"""
function singularF1(a::FT, b::FT, c::FT, d::FT) where{FT<:AbstractFloat}
    s   = (a + b + c + d)/2
    return  -4* (   1/a*log(1 - a/s) + 
                    1/b*log(1 - b/s) + 
                    1/c*log(1 - c/s) +
                    1/c*log(1 - d/s)  )/3
end

@doc raw"""
    singularF21(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}

计算边长为`a, b, c`，面积平方为`area2`的三角形重合时的奇异性F2项，即
``\int{\int{\frac{\boldsymbol{\rho}_{m}\cdot\boldsymbol{\rho}_{n}}{R}}dS}``
的解析值，该函数处理 m==n 即基函数自作用的情况。
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

@doc raw"""
    singularF21(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}

计算边长为`a, b, c`，面积平方为`area2`的三角形重合时的奇异性F2项，即
``\int{\int{\frac{\boldsymbol{\rho}_{m}\cdot\boldsymbol{\rho}_{n}}{R}}dS}``
的解析值，该函数处理 m!=n 即同一三角形的不同基函数作用的情况。
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

@doc raw"""
    faceSingularityIg(rgt::AbstractVector{FT}, tris::TriangleInfo{IT, FT}, area::FT, facen̂::AbstractVector{FT}) where {IT<:Integer, FT<:Real}

计算场点`rgt`在源三角形`tris`上的奇异性，`tris`的面积为`area`，外法向量为`facen̂`。
计算结果为：
```math
\begin{aligned}
I_{gS}  &= \int{g(R)dS'}\\
        &= \sum_{n=0}^{SglrOrder}{coeffgreen(n)I^{n-1}_{RS}}\\
I^{n}_{RS}  &= \int{R^{n}dS'}
\end{aligned}
```
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

@doc raw"""
    faceSingularityIg(rgt::AbstractVector{FT}, tris::Tris4Tetra{IT, FT}, area::FT, facen̂::AbstractVector{FT}) where {IT<:Integer, FT<:Real}

计算场点`rgt`在源三角形`tris`（该三角形为组成四面体的某一面）上的奇异性，`tris`的面积为`area`，外法向量为`facen̂`。
计算结果为：
```math
\begin{aligned}
I_{gS}  &= \int{g(R)dS'}\\
        &= \sum_{n=0}^{SglrOrder}{coeffgreen(n)I^{n-1}_{RS}}\\
I^{n}_{RS}  &= \int{R^{n}dS'}
\end{aligned}
```
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

@doc raw"""
    faceSingularityIg(rgt::AbstractVector{FT}, polys::ST, area::FT, 
    facen̂::AbstractVector{FT}) where {IT<:Integer, FT<:Real, ST<:SurfaceCellType{IT, FT}}

计算场点`rgt`在多边形`polys`上的奇异性，`polys`的面积为`area`，外法向量为`facen̂`。
计算结果为：
```math
\begin{aligned}
I_{gS}  &= \int{g(R)dS'}\\
        &= \sum_{n=0}^{SglrOrder}{coeffgreen(n)I^{n-1}_{RS}}\\
I^{n}_{RS}  &= \int{R^{n}dS'}
\end{aligned}
```
"""
function faceSingularityIg(rgt::AbstractVector{FT}, polys::ST, area::FT, 
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
    # dts     =   facen̂ ⋅ (rgt .- view(polys.vertices, :, 1))
    dts  =   zero(FT)
    for ii in 1:3
        dts += facen̂[ii] * (rgt[ii] - polys.vertices[ii, 1])
    end
    # 第 gi 个投影点
    r0gi    =   rgt .- dts .* facen̂
    
    # 距离的绝对值、平方
    dtsAbs  =   abs(dts)
    dts²    =   dtsAbs^2
    # 边数
    edgeN   =   length(polys.edgel)
    
    # 预分配临时变量以加速
    p02jvec =   zero(MVec3D{FT})

    # 对三个边循环
    for edgej in 1:edgeN
        # 构成该边的第一个点
        edgeNodei⁻  =    polys.vertices[:, edgej]
        # 该边边长
        lⱼ          =   abs(polys.edgel[edgej])
        # 每个边的局部坐标下的 lⱼ⁺ \lⱼ⁻
        # lⱼ⁻     =   (edgeNodei⁻ .- r0gi) ⋅ polys.edgev̂[:, edgej]
        lⱼ⁻  =   zero(FT)
        for ii in 1:3
            lⱼ⁻ += (edgeNodei⁻[ii] - r0gi[ii]) * polys.edgev̂[ii, edgej]
        end
        lⱼ⁺     =   lⱼ⁻ + lⱼ
        # 投影点 r0gi 到各个边的垂足的向量 p02jvec
        # p02jvec =   edgeNodei⁻ .- lⱼ⁻ * polys.edgev̂[:, edgej] .- r0gi
        for ii in 1:3
            p02jvec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * polys.edgev̂[ii, edgej] - r0gi[ii]
        end
        # 该向量长的平方、该向量长
        p02jl   =   p02jvec ⋅ polys.edgen̂[:, edgej]
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

@doc raw"""
    faceSingularityIgIvecg(rgt::AbstractVector{FT}, polys::ST, area, 
        facen̂::AbstractVector) where {FT<:Real, ST<:SurfaceCellType{IT, FT}}

计算场点`rgt`在多边形`polys`上的奇异性，`polys`的面积为`area`，外法向量为`facen̂`。
计算结果为：
```math
\begin{aligned}
I_{gS}  &= \int{g(R)dS'}\\
        &= \sum_{n=0}^{SglrOrder}{coeffgreen(n)I^{n-1}_{RS}}\\
\boldsymbol{I}_{gS}  &= \int{\boldsymbol{R}g(R)dS'}\\
        &= \sum_{l_j}{\hat{\bm{u}}_j \sum_{n=0}^{SglrOrder}{\frac{coeffgreen(n)}{n+1}I^{n-1}_{lr}}} + d\bm{n}I^{n-1}_{gS}\\
\end{aligned}
```
"""
function faceSingularityIgIvecg(rgt::AbstractVector{FT}, polys::ST, area, 
    facen̂::AbstractVector) where {FT<:Real, ST<:SurfaceCellType}

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
    # dts     =   facen̂ ⋅ (rgt .- view(polys.vertices, :, 1))
    dts  =   zero(FT)
    for ii in 1:3
        dts += facen̂[ii] * (rgt[ii] - polys.vertices[ii, 1])
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
        edgeNodei⁻  =    polys.vertices[:, MoM_Basics.EDGEVmINTriVsID[edgej]]
        # 该边边长
        lⱼ          =   abs(polys.edgel[edgej])
        # 每个边的局部坐标下的 lⱼ⁺ \lⱼ⁻
        # lⱼ⁻     =   (edgeNodei⁻ .- r0gi) ⋅ polys.edgev̂[:, edgej]
        lⱼ⁻  =   zero(FT)
        for ii in 1:3
            lⱼ⁻ += (edgeNodei⁻[ii] - r0gi[ii]) * polys.edgev̂[ii, edgej]
        end
        lⱼ⁺     =   lⱼ⁻ + lⱼ
        # 投影点 r0gi 到各个边的垂足的向量 p02jvec
        # p02jvec =   edgeNodei⁻ .- lⱼ⁻ * polys.edgev̂[:, edgej] .- r0gi
        for ii in 1:3
            p02jvec[ii] = edgeNodei⁻[ii] - lⱼ⁻ * polys.edgev̂[ii, edgej] - r0gi[ii]
        end
        # 该向量长的平方、该向量长
        p02jl   =   p02jvec ⋅ polys.edgen̂[:, edgej]
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
            @views IvecSg   .+=  polys.edgen̂[:, edgej] .* Cdvnp1I
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