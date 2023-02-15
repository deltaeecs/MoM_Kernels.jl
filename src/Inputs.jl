"""
用于输入频率参数，修改其它仿真参数的函数
输入积分方程类型参数
"""
function inputParameters(;args...)
    inputBasicParameters(;args...)
    modiSingularityRelatedConsts!()
    return nothing
end
