# 快速算法参数

## 叶层正方体盒子边长:leaf_cube_side_length
# const LEAFCUBESIZE  =   0.23*Params.λ_0

"""
创建可变参数类型以在频率更改时对应更改 MLFMA 的相关参数
"""
mutable struct MLFMAParamsType{FT<:AbstractFloat}
    LEAFCUBESIZE ::FT
end

const MLFMAParams   =   MLFMAParamsType{typeof(Params.frequency)}(0.23*Params.λ_0)

function modiMLFMAParams!(λ_0::FT) where {FT<:AbstractFloat}
    MLFMAParams.LEAFCUBESIZE = 0.23*λ_0
end

# ## 严格使用此长度（:strict）还是根据目标自适应（：recomend能适当修改叶层边长使得层数更合适 ）？
# const LEAFCUBESIZEMODE  =   :recomend
# MLFMA精度
const NBDIGITS      =   IntDtype(3)


"""
拉格朗日局部插值每个方向上的插值点个数, 要设置为偶数
"""
const NLocalInterpolation = 6
(isodd(NLocalInterpolation) | NLocalInterpolation <= 0) && throw(ArgumentError("NLocalInterpolation 局部插值点数应为偶数"))

