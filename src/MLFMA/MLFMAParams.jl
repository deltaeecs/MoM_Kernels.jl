# 快速算法参数

## 叶层正方体盒子边长:leaf_cube_side_length
# const LEAFCUBESIZE  =   0.23*Params.λ_0

"""
创建可变参数类型以在频率更改时对应更改 MLFMA 的相关参数
"""
Base.@kwdef mutable struct MLFMAParamsType{FT<:AbstractFloat}
    NBDIGITS::Int = 3# MLFMA精度
    LEAFCUBESIZE ::FT = 0.23*Params.λ_0
    InterpolationMethod::Symbol = :Lagrange2Step
end

const MLFMAParams   =   MLFMAParamsType{typeof(Params.frequency)}(LEAFCUBESIZE = 0.23*Params.λ_0)

function modiMLFMAParams!(λ_0::FT) where {FT<:AbstractFloat}
    MLFMAParams.LEAFCUBESIZE = 0.23*λ_0
end

function set_Interpolation_Method!(method)
    MLFMAParams.InterpolationMethod = method
    nothing
end

function get_Interpolation_Method(method)
    if method == :Lagrange2Step
        return LagrangeInterpInfo
    elseif method == :Lagrange1Step
        return LagrangeInterp1StepInfo
    else
        throw("插值方法选择出错")
    end
end