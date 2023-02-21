"""
从 IncompleteLU.jl[https://github.com/haampie/IncompleteLU.jl.git] 
包实现ilu, 再次封装是因为要加入一些判断
"""
function iluPrecondition(A, level; τ = 1e-3)

    if typeof(A) <: Transpose
        @warn "近场矩阵采用了转置以加快计算速度，不再适用 ilu 算法，改回 SAI!"
        return sparseApproximateInversePl(A, level)
    else
        @clock "计算ILU" begin
            re = ilu(A; τ = τ)
        end
        return re
    end

end