"""
    inputParameters(;args...)

用于输入仿真参数，并修改奇异性处理中频率相关常量。
详见 [`MoM_Basics.inputBasicParameters`](@ref) 和 [`modiSingularityRelatedConsts!`](@ref)。
"""
function inputParameters(;args...)
    inputBasicParameters(;args...)
    set_leafCubeSize!()
    modiSingularityRelatedConsts!()
    return nothing
end
