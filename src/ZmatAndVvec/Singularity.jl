
# 计算奇异性的阶数 ( 即格林函数泰勒展开的阶数，越高阶越精确 )
const SglrOrder =   15

@doc raw"""
    greenfunc_star(pa::Vec3D{T}, pb::Vec3D{T}[; k=Params.K_0, taylorOrder = SglrOrder]) where {T<:AbstractFloat}

归一化格林函数 ``4πG(R)`` 泰勒展开后去奇异项：
```math
g^{*}(R)    = \frac{e^{-jkR}}{R} - \frac{1}{R}
            = \sum_{n=1}^{SglrOrder}\frac{{-jk}^{n}}{n!}R^{n-1}.
```
"""
function greenfunc_star(ra::Vec3D{T}, rb::Vec3D{T}; k=Params.K_0, taylorOrder = SglrOrder) where {T<:AbstractFloat}
    # 两点距离
    R   =   dist(ra, rb)
    return greenfunc_star(R; k = k, taylorOrder=taylorOrder)
end
function greenfunc_star(R::T; k=Params.K_0, taylorOrder = SglrOrder) where {T<:AbstractFloat}
    # 一些循环中临时变量
    minusJk =   -1im*k
    # 先初始化为1阶解用于累加
    g_star  =   minusJk
    temp0   =   minusJk*R
    temp1   =   minusJk
    
    # 格林函数计算
    taylorOrder >= 2 && @inbounds for i in 2:taylorOrder
        temp1   *=  temp0/i
        g_star  +=  temp1
    end #for

    return g_star
end

@doc raw"""
    coeffgreen(n::Integer)

归一化格林函数 (不包括``\frac{1}{4π}``项) 的展开系数函数，从 0 阶 到 n 阶：
```math
\begin{aligned}
g(R)    &= \frac{e^{-jkR}}{R} = \sum_{n=0}^{SglrOrder}coeffgreen(n)R^{n-1}\\
coeffgreen(n)  &=   \frac{{-jk}^{n}}{n!}\\
\end{aligned}
```
"""
coeffgreen(n::Integer) = (-Params.JK_0)^n/factorial(n)

## 用到的常量
"""
计算面奇异性时用的两个系数。
"""
const SSCg          =   OffsetVector(MVector{SglrOrder+1, Precision.CT}(map(n -> coeffgreen(n), 0:SglrOrder)), -1)
const SSCgdivnp1    =   OffsetVector(MVector{SglrOrder+1, Precision.CT}(map(n -> coeffgreen(n)/(n+1), 0:SglrOrder)), -1)

"""
计算体奇异性时用到的一些系数。
"""
const SSCgdivnp2    =   OffsetVector(MVector{SglrOrder+1, Precision.CT}(map(n -> coeffgreen(n)/(n+2), 0:SglrOrder)), -1)
const VSC₁ⁿ   =   OffsetArray(zero(MVector{SglrOrder+1, Precision.CT}), -2)
const VSC₂ⁿ   =   OffsetArray(zero(MVector{SglrOrder-1, Precision.CT}), -2)

"""
这一项面、体奇异性都用到了，处理对格林函数梯度求积时的奇异性时用到的。
"""
const VSC₃ⁿ   =   OffsetArray(zero(MVector{SglrOrder+1, Precision.CT}), -2)

"""
    setVSC₁₂₃ⁿ!()

计算体奇异性三个系数。
"""
function setVSC₁₂₃ⁿ!()
    for idx in eachindex(VSC₁ⁿ)
        VSC₁ⁿ[idx] =  coeffgreen(idx+1)/(idx + 3)
        VSC₃ⁿ[idx] =  coeffgreen(idx+1)
    end
    for idx in eachindex(VSC₂ⁿ)
        VSC₂ⁿ[idx] =  Params.JK_0*coeffgreen(idx+2) + coeffgreen(idx+3)
    end
    nothing
end

# 导入面、体奇异性处理函数
include("Singularity/FaceSingularity.jl")
include("Singularity/VolumeSingularity.jl")


"""
用于输入参数（特指频率）改变时的更改相关常数项
"""
function modiSingularityRelatedConsts!()
    SSCg          .=   OffsetVector(MVector{SglrOrder+1}(map(n -> coeffgreen(n), 0:SglrOrder)), -1)
    SSCgdivnp1    .=   OffsetVector(MVector{SglrOrder+1}(map(n -> coeffgreen(n)/(n+1), 0:SglrOrder)), -1)
    SSCgdivnp2    .=   OffsetVector(MVector{SglrOrder+1}(map(n -> coeffgreen(n)/(n+2), 0:SglrOrder)), -1)
    setVSC₁₂₃ⁿ!()
    nothing
end
