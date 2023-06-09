# SAI
include("SAI.jl")
include("SAIChunks.jl")
# ilu
include("ILU.jl")

using IncompleteLU:ILUFactorization

@doc """
    size(operator::T) where {T<:ILUFactorization}

提供 ilu 的算子 size 函数。
"""
size(operator::T) where {T<:ILUFactorization} = (size(operator.L, 1), size(operator.U, 2))
function size(operator::T, idx::I) where {T<:ILUFactorization, I<:Integer}
    idx > 3 && throw(ArgumentError("索引错误，ilu 算子为 2 维的。"))
    size(getfield(operator, idx), idx)
end

"""
    eltype(opt::ILUFactorization)

提供 ilu 的算子 eltype 函数。
"""
eltype(opt::ILUFactorization) = eltype(opt.L)
