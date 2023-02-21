"""
展示
"""
function record_memorys(Zopt)

    memory["近场矩阵"]  = Base.summarysize(Zopt.Znear)
    memory["八叉树"]    = Base.summarysize(Zopt.octree)
    memory["网格信息"]   = Base.summarysize(Zopt.vsCellsInfo)
    memory["叶层聚合项"] = Base.summarysize(Zopt.aggSBF)
    
    nothing
end


"""
    restore_infos()
    记录各部分内存和各阶段计算时间。
TBW
"""
function restore_infos(;targetfile = joinpath(SimulationParams.resultDir, "InputArgs.txt"), mode = "a+")

    open(targetfile, mode)  do f
        show_memory_time(f)
    end

end