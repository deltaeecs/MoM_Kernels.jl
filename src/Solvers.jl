## 本文件用于构建求解器
"""
迭代求解器选择
"""
function iterSolverSet(solverType::Symbol)::Function
    getproperty(IterativeSolvers, solverType)
end

"""
保存电流系数
"""
function saveCurrent(ICurrent; str = "")
    jldsave(SimulationParams.resultDir*"ICurrent$str.jld2"; ICurrent = ICurrent)
end

"""
读取电流系数
"""
function loadCurrent(filename)
    occursin("jld2", filename) ? load("$filename", "ICurrent") : load("$filename.jld2", "ICurrent")
end


"""
矩阵方程
Ax=b
复合求解函数
输入值：
A::Matrix{T}, b::Vector{T}
solverT::Symbol  求解器类型
"""
function solve(A::LinearMapType{T}, b::AbstractVector{T};
    solverT::Symbol=:gmres, Pl = Identity(), Pr = Identity(), rtol = 1e-3,
    maxiter = 1000, str = "", restart = 200) where{T<:Number}
    # 直接求解
    solverT  == :direct  && begin
        println("Solving matrix function with LUD.")
        try
            return (A\b, nothing)
        catch
            println("查看输入是否为矩阵！")
            return nothing, nothing
        end
    end  

    FT = real(T)
    # 迭代求解器
    solver      =   iterSolverSet(solverT)
    # 残差阈值
    resnorm0    =   FT(norm( Pl \ b))
    resnormtol  =   FT(rtol*resnorm0)

    # 迭代求解
    println("Solving matrix function with iterate solver $solverT, \nwith initial resnorm $resnorm0.")
    ICoeff, ch       =   solver(A, b; restart = restart, abstol = resnormtol, Pl = Pl, Pr = Pr,  log = true, verbose=true, maxiter = maxiter)
    saveCurrent(ICoeff; str = str)

    return ICoeff, ch
end

"""
矩阵方程
Ax=b
复合求解函数
输入值：
A::Matrix{T}, b::Vector{T}
solverT::Symbol  求解器类型
"""
function solve!(A::LinearMapType{T}, x::AbstractVector{T}, b::AbstractVector{T}; 
    solverT::Symbol = :gmres!, Pl = Identity(), Pr = Identity(), rtol = 1e-3, 
    maxiter = 1000, str = "", restart = 200) where{T<:Number}
    # 直接求解
    solverT  == :direct  && begin
        println("Solving matrix function with LUD.")
        try
            return (A\b, nothing)
        catch
            println("查看输入是否为矩阵！")
            return nothing, nothing
        end
    end  

    FT = real(T)
    # 迭代求解器
    solver      =   iterSolverSet(solverT)
    # 残差阈值
    resnorm0    =   FT(norm( Pl \ b))
    resnormtol  =   FT(rtol*resnorm0)

    # 迭代求解
    println("Solving matrix function with iterate solver $solverT, with initial resnorm $resnorm0. ")
    # try 
    ICoeff, ch       =   solver(x, A, b; restart = restart, abstol = resnormtol, Pl = Pl, Pr = Pr,  log = true, verbose=true, maxiter = maxiter)
    saveCurrent(ICoeff; str = str)

    # 相对残差结果
    relresnorm  =   ch.data[:resnorm] / resnorm0

    # 命令行绘图
    SimulationParams.SHOWIMAGE  &&  convergencePlot(relresnorm)

    # 将相对残差写入文件
    open(joinpath(SimulationParams.resultDir, "$(solverT)_ch$str.txt"), "w") do io
        for resi in relresnorm
            write(io, "$resi\n" )
        end
    end

    return ICoeff, ch

end

"""
计算完成后绘制收敛曲线
"""
function convergencePlot(resnorm::Vector{FT}) where{FT<:Real}

    figCnvergence = lineplot(resnorm, yscale = :log10, ylim = (minimum(resnorm), 1), xlabel = "Epoch", title = "Relative ResNorm - Epoch")

    SimulationParams.SHOWIMAGE  &&  display(figCnvergence)

    return 

end