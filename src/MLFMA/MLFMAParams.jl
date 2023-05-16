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
    useCSR::Bool = false
end

"""
多层快速多极子的可调参数
"""
const MLFMAParams   =   MLFMAParamsType{typeof(Params.frequency)}(LEAFCUBESIZE = 0.23*Params.λ_0)

"""
设置叶层盒子边长
"""
function set_leafCubeSize!(cubel::FT = 0.25Params.λ_0) where {FT<:AbstractFloat}
    MLFMAParams.LEAFCUBESIZE = cubel
end

"""
获取叶层盒子边长
"""
function get_leafCubeSize()
    MLFMAParams.LEAFCUBESIZE
end

"""
设置插值算法
"""
function set_Interpolation_Method!(method)
    MLFMAParams.InterpolationMethod = method
    nothing
end

@doc """
    use_CSR()
    use_CSC()
设置近场阻抗矩阵是否采用CSR（CSC转置）。
"""
use_CSR() = setfield(MLFMAParams, :useCSR, true)
use_CSC() = setfield(MLFMAParams, :useCSR, false)

@doc """
    get_Interpolation_Method(method::Symbol)
    get_Interpolation_Method(method::Union{Val{:Lagrange2Step}, Val{:Lagrange1Step}})

获取插值算法。
"""
get_Interpolation_Method(method::Symbol) = get_Interpolation_Method(Val(method))
function get_Interpolation_Method(method::Union{Val{:Lagrange2Step}, Val{:Lagrange1Step}})
    if method == Val(:Lagrange2Step)
        return LagrangeInterpInfo
    elseif method == Val(:Lagrange1Step)
        return LagrangeInterp1StepInfo
    else
        throw("插值方法选择出错")
    end
end