"""
在叶层从基函数向盒子聚合
"""
function aggOnTargetLevel(level, aggSBF)
    # 叶层盒子
    cubes   =   level.cubes
    # 叶层聚合项
    aggS    =   copy(level.aggS)
    aggS   .=   0
    # 对盒子循环计算
    @threads for iCube in eachindex(cubes)
        # 盒子信息
        cube    =   cubes[iCube]
        # 基函数区间
        bfInterval  =   cube.bfInterval
        # 往盒子中心聚合
        @inbounds for n in bfInterval
            @views aggS[:, :, iCube]   .+=  aggSBF[:, :, n]
        end #n
    end #iCube
    return aggS
end

levels = octree.levels

for level in values(levels)
    setGeoIDsInLeafCubes!(level, rwgsInfo)
end

levels = octree.levels
targetLevel =    levels[2];


ICoeffone = ones(Complex{Float64}, length(V))
calZI*ICoeffone


# 目标聚合向
aggTargetBF, _    =   aggSBFOnLevel(targetLevel, trianglesInfo, rwgsInfo)

aggTarget = aggOnTargetLevel(targetLevel, aggTargetBF)

aggInterped = targetLevel.aggS

εs = zeros(size(aggInterped))
for i in 1:length(targetLevel.cubes)
    εs[:, 1, i] .= abs.(aggInterped[:, 1, i] .- aggTarget[:, 1, i]) ./ maximum(abs.(aggTarget[:, 1, i]))
    εs[:, 2, i] .= abs.(aggInterped[:, 2, i] .- aggTarget[:, 2, i]) ./ maximum(abs.(aggTarget[:, 2, i]))
end

εsmean = mean(εs, dims = [1, 3])

Plots.histogram(εsmean[:])

# 5-2
# 3.5446722441350035e-5
# 4.06783e-5  3.02151e-5

# 5-3
# 5.7474235314763334e-5
# 6.67667e-5  4.81817e-5

# 5-4
# 4.380214653525037e-5
# 4.96995e-5  3.79048e-5 -#