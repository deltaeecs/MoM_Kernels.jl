module MoM_Kernels

# 依赖于 MoM_Basics 提供的各种类型
using MoM_Basics

# 其他依赖
using ProgressMeter, JLD2
using StaticArrays, OffsetArrays, SparseArrays
using .Threads, ThreadsX, FLoops, FoldsThreads
using LinearAlgebra, Statistics
using FastGaussQuadrature, SpecialFunctions, GSL, LegendrePolynomials
using IncompleteLU, IterativeSolvers


export  getOctreeAndReOrderBFs!,
        calZnearCSC, getImpedanceMatrix,
        getExcitationVector, getImpedanceMatAndExciteV,
        MLMFAIterator,
        SAIPrec, SAIChunkPrec,
        sparseApproximateInversePl, sparseApproximateInversePr,
        solve, solve!


# MLFMA参数
include("MLFMA/MLFMAParams.jl")

# 采用直接算法或者快速算法
include("ZmatAndVvec/ZmatVvec.jl")
include("DirectAlgorithm.jl")
include("FastAlgorithm.jl")

# 求解函数
include("Solvers.jl")

end
