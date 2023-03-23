## 本文件存储的函数用于计算八叉树上的积分要用到的一些量，如截断项数、角谱空间采样点、采样权重等

"""
该函数计算八叉树各层截断项数
输入为本层最小盒子的边长
"""
function truncationLCal(cubel::FT) where {FT<:Real}
    L = floor(Int, cubel*real(Params.K_0)*sqrt(3) + 2.16*MLFMAParams.NBDIGITS^(2.0/3.0)*(cubel*real(Params.K_0))^(1/3))
    return L
end

## 球面多极子（采样点）信息抽象类型
abstract type PolesInfo{FT<:AbstractFloat} end

## 插值信息抽象类型
abstract type InterpInfo{IT<:Integer, FT<:Real} end

## 拉格朗日插值
include("LagrangeInterpolation.jl")